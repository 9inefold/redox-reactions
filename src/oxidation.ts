import { type ElementName, type RawElement, type Element, type Molecule, type MoleculePart, hashMolecule } from "./peg";
import { MoleculeData } from "./lewis";
import { periodicTable } from './elements';
import { type LaxRecord } from "./object-utils";

export type OxidationStates = LaxRecord<ElementName, number>;

// Known oxidation states for common elements/ions
const FIXED_OXIDATION_STATES: Record<string, number> = Object.freeze({
  F:  -1,
  Li: +1, Na: +1, K:  +1, Rb: +1, Cs: +1, // Group 1
  Ca: +2, Mg: +2, Ba: +2, Sr: +2, Be: +2, // Group 2
  Al: +3, // Group 13
  Ag: +1,
  Zn: +2,
});

const KNOWN_STATES: Record<string, OxidationStates> = {
  'H1O1{-1*0': { H: +1, O: -2 },
  'H2O1{0*0': { H: +1, O: -2 },
  'H2O2{0*0': { H: +1, O: -1 },
  'O4S1{-2*0': { O: -2, S: +6 },
  'O3S1{-2*0': { O: -2, S: +4 },
  'C1O3{-2*0': { C: +4, O: -2 },
  'N1O3{-1*0': { N: +5, O: -2 },
  'N1O2{-1*0': { N: +3, O: -2 },
  'Cl1O4{-1*0': { Cl: +7, O: -2 },
  'Cl1O3{-1*0': { Cl: +5, O: -2 },
  'Cl1O2{-1*0': { Cl: +3, O: -2 },
  'Cl1O1{-1*0': { Cl: +1, O: -2 },
  'C2H3O2{-1*0': { C: +3, H: +1, O: -2 }
};

function checkElement(element: Element) {
  if (element.name in periodicTable)
    return;
  throw new Error(`Unknown element: ${element.name}`);
}

export function flattenElements(data: (Element | MoleculePart)[]): RawElement[] {
  const result: RawElement[] = [];
  for (const part of data) {
    if (part.element) {
      checkElement(part);
      result.push({
        count: part.count,
        name: part.name as ElementName
      });
    } else {
      for (const el of part.data) {
        checkElement(el);
        result.push({
          count: el.count * part.count,
          name: el.name as ElementName
        });
      }
    }
  }
  return result;
}

export function getOxidationStates(molecule: Molecule): OxidationStates | undefined {
  const data = new MoleculeData(molecule);
  const result: OxidationStates = {};

  if (data.unique === 1) {
    // Monoatomic
    const A = data.central;
    result[A] = molecule.charge / data.counts[A]!;
    return result;
  }

  if (data.id in KNOWN_STATES)
    return KNOWN_STATES[data.id];

  const elts = Object.keys(data.counts) as ElementName[];
  const knowns: ElementName[] = [];
  const unknowns: ElementName[] = [];

  for (const el of elts) {
    if (el in FIXED_OXIDATION_STATES)
      knowns.push(el);
    else
      unknowns.push(el); 
  }

  function popOutOH(OH: ElementName) {
    const OHidx = unknowns.findIndex(unk => unk === OH);
    unknowns[OHidx] = unknowns[unknowns.length - 1];
    unknowns.pop();
    //knowns.push(OH);
  };

  const uniqueMetals = elts.reduce((prev, str) => {
    return periodicTable[str].metal ? prev + 1 : prev;
  }, 0);

  const hasH = 'H' in data.counts;
  const hasO = 'O' in data.counts;

  function calculateOHOx(): number {
    if (hasH && !hasO) {
      const hasB = ('B' in data.counts) ? 1 : 0;
      const metals = uniqueMetals + hasB;
      // Too complex!
      if (metals > 1) return 0;
      else if (metals >= 1 && elts.length > 2) return 0;
      const ox = metals === 0 ? +1 : -1;
      popOutOH('H');
      result['H'] = ox;
      return ox * data.counts['H']!;
    } else if (hasO && !hasH) {
      popOutOH('O');
      result['O'] = -2;
      return -2 * data.counts['O']!;
    }
    // Nothing extra
    return 0;
  };

  // Calculate just O/H
  const OHCharge = calculateOHOx();
  // Calculate the rest
  const knownCharge = knowns.reduce((prev, str) => {
    const ox = FIXED_OXIDATION_STATES[str]!;
    result[str] = ox;
    return prev + (ox * data.counts[str]!);
  }, 0) + OHCharge;

  // Check for the simplest case
  if (unknowns.length === 0 && knownCharge === molecule.charge) {
    return result;
  }
  // Check for a simple algabraic solution
  if (unknowns.length === 1) {
    // Solve x + K = C
    const unk = unknowns[0];
    const x = molecule.charge - knownCharge;
    // Add in unknown
    result[unk] = x / data.counts[unk]!;
    return result;
  }

  if (uniqueMetals > 2)
    throw Error('Too many metals!');

  // Checks for things like H2(SO4)
  function tryStrippingEl(name: ElementName): OxidationStates | undefined {
    const max = periodicTable[name].numberOfValence;
    if (!max) return undefined;
    // Try to match against a known anion
    let removedH = 0;
    const stripped = molecule.data.filter(e => {
      if (!e.element || e.name !== name) return true;
      // Count the number removed
      removedH += e.count;
      return false;
    });
    // Ensure hydrogen actually got stripped
    if (removedH === 0)
      return undefined;
    // Generate a new hash for stripped molecule
    let ox = undefined;
    let valence = 1;
    for (; valence <= max; ++valence) {
      const hash = hashMolecule({
        charge: molecule.charge - (removedH * valence),
        hydrate: molecule.hydrate,
        data: stripped
      });
      ox = KNOWN_STATES[hash];
      if (ox !== undefined)
        break;
    }
    // Check if our new one exists
    if (ox === undefined)
      return undefined;
    // Check if it already has H
    if (name === 'H' && 'H' in ox) {
      // Ensure theres no weirdness here
      if (ox['H'] !== +1)
        return undefined;
      return { ...ox };
    }
    // No hydrogen, add ours
    let out = { ...ox };
    out[name] = valence;
    return out;
    //return { H: removedH, ...ox };
  };

  if (hasH && hasO) {
    const alternate = tryStrippingEl('H');
    if (alternate !== undefined)
      return alternate;
  }

  {
    const metal = unknowns.find(m => periodicTable[m].metal);
    if (metal !== undefined) {
      const alternate = tryStrippingEl(metal);
      if (alternate !== undefined)
        return alternate;
    }
  }

  //data.log('NEED LEWIS STRUCTURE!!!');
  return undefined;
}
