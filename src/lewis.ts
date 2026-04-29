import type { RawElement, Element, MoleculePart, Molecule, ElementName } from "./peg";
import { hashMolecule } from "./peg";
import { flattenElements } from './oxidation';
import { periodicTable } from "./elements";

export class MoleculeData {
  readonly molecule: Molecule;
  readonly flat: RawElement[];
  readonly counts: ElementRecord<number>;
  readonly unique: number;
  private _id?: string;

  private _central?: ElementName;
  private _multiCenter?: boolean;
  private _totalVE: number = 0;
  private _totalBE: number = 0;
  private _totalE: number = 0;

  private _warnings: string[] = [];

  constructor(molecule: Molecule) {
    this.molecule = molecule;
    this.flat = flattenElements(molecule.data);
    this.counts = collectElementCounts(this.flat);
    this.unique = Object.keys(this.counts).length;

    this._setupEs();
  }

  log(...args: any[]) {
    if (args.length !== 0)
      console.log(`${this.molecule.text}: `, ...args);
    else
      console.log();
  }

  warn(message: string) {
    this._warnings.push(message);
  }

  get warnings(): string[] {
    return this._warnings;
  }

  get id(): string {
    if (this._id === undefined)
      this._id = hashMolecule(this.molecule);
    return this._id;
  }

  get central(): ElementName {
    if (this._central === undefined)
      this._setupCentral();
    return this._central!;
  }

  get multiCenter(): boolean {
    if (this._multiCenter === undefined)
      this._setupCentral();
    return this._multiCenter!;
  }

  get totalVE(): number {
    if (this._totalVE !== undefined)
      return this._totalVE;
    // Calculate totalVE
    this._totalVE = 0;
    for (const [name, count] of Object.entries(this.counts))
      this._totalVE += getValence(name as ElementName) * count;
    // add electrons for negative charge, remove for positive
    this._totalVE -= this.molecule.charge;
    return this._totalVE;
  }

  private _setupCentral() {
    this._central = pickCentralAtom(this.flat);
    this._multiCenter = this.counts[this._central]! > 1;
  }

  private _setupEs() {
    // Calculate totalVE/totalBE
    for (const [name, count] of Object.entries(this.counts)) {
      this._totalVE += getValence(name as ElementName) * count;
      this._totalE += name === 'H' ? 2 : 8;
    }
    // add electrons for negative charge, remove for positive
    this._totalVE -= this.molecule.charge;
    this._totalBE = this._totalVE - this._totalE;
  }
};

/** 1 = single, 2 = double, 3 = triple */
export type BondOrder = 1 | 2 | 3;

/** A single atom node in the Lewis graph */
export interface Atom {
  id: number;
  name: ElementName;
  lonePairs: number;
  /** Formal charge = V - E - B/2 */
  formalCharge: number;
  /** How many bonds this atom participates in (counting multiplicity) */
  bondOrder: number;
};

/** A bond edge between two atoms */
export interface Bond {
  id: number;
  from: number;   // atom id
  to: number;     // atom id
  order: BondOrder;
};

/** The full Lewis structure */
export interface LewisStructure {
  formula: string;
  charge: number;            // molecular charge (0 = neutral)
  totalValenceElectrons: number;
  atoms: Atom[];
  bonds: Bond[];
  /** True if all atoms satisfy octet/duet rule with minimal formal charges */
  isValid: boolean;
  formalChargeSum: number;
  description: string;
  warnings: string[];
};

type ElementRecord<T> = Partial<Record<ElementName, T>>;

const MONOCARBON: Element = Object.freeze({ element: true, count: 1, name: 'C' });

function getValence(name: ElementName): number {
  return periodicTable[name].numberOfValence ?? 0;
}

////////////////////////////////////////////////////////////////////////////////
// Kekule structure
////////////////////////////////////////////////////////////////////////////////

type KekuleElement = RawElement & { valence: number };

type KekulePart = (Element | MoleculePart) & { valence?: number };

/** Valid groups */
export type KekuleGroup = KekulePart[];

/**
 * Determine if the given formula is a kekule structure.
 * This will give us a hint as to the structure of the molecule.
 */
function getKekuleStructure(molecule: Molecule): KekuleGroup[] | undefined {
  let structures: KekuleGroup[] = [];
  let currCount = 0;
  let current: KekuleGroup = [];
  let isFirst = true, isEnd = false;

  function pushCurr() {
    if (currCount > 3) {
      if (isFirst)
        isFirst = false;
      else
        isEnd = true;
    }
    currCount = 0;
    structures.push(current);
  };

  function fail(message: string): undefined {
    console.log(`${molecule.text}: not kekule, ${message}`);
    return undefined;
  }

  for (const el of molecule.data) {
    if (isEnd) return fail('bonds after cap bond');
    // Max 3 bonds on the central path
    if (currCount > 4) return fail(`carbon overbonded (${currCount})`);
    const count = current.length;

    // Handle the first element
    if (isFirst && count === 0) {
      if (el.element) {
        // Ensure this is a valid first bond
        if (el.name !== 'C' || el.count > 1)
          return fail(`expected C1 cap, got ${el.name}${el.count}`);
        current.push(el);
        ++currCount;
      } else {
        const subEl = el.data[0];
        // Ensure this subbond is valid
        if (subEl.name !== 'C' || subEl.count > 1)
          return fail(`expected C1 cap, got ${subEl.name}${subEl.count}`);
        // Add whole group
        current.push(...el.data);
        currCount = el.data.reduce((sum, val) => sum + val.count, 0);
      }
      continue;
    }

    // Handle new start
    if (el.element && el.name === 'C') {
      pushCurr();
      for (let i = 0; i < el.count - 1; ++i)
        structures.push([MONOCARBON]);
      current.push(MONOCARBON);
      ++currCount;
      continue;
    }

    // Handle other cases
    current.push(el);
    currCount += el.count;
  }

  if (currCount > 4)
    // Fumbled at the finish line
    return fail(`cap carbon overbonded (${currCount})`);
  
  pushCurr();
  return structures;
}

const VALENCE_OPTIONS: ElementRecord<number[]> = {
  H:  [1],
  C:  [4],
  N:  [3, 5],
  O:  [2],
  P:  [3, 5],
  S:  [2, 4, 6],
  F:  [1],
  Cl: [1, 3, 5, 7],
  Br: [1, 3, 5],
  I:  [1, 3, 5, 7],
  B:  [3],
  Si: [4],
  Se: [2, 4, 6],
  As: [3, 5],
};

/**
 * Returns the list of valence options for an element name.
 * Falls back to the periodic table's `valence` field if present,
 * or [0] as a last resort.
 */
function getValenceOptions(name: ElementName): number[] {
  if (VALENCE_OPTIONS[name]) return VALENCE_OPTIONS[name]!;
  return [getValence(name)]
}

/**
 * Charge adjustment: a molecule with charge q has q fewer electrons than
 * neutral, shifting the effective valence sum by -q (cations lose electrons,
 * anions gain them).  We redistribute the charge across atoms that have
 * multiple valence options by trying different valence assignments.
 *
 * For simplicity we adjust the atom with the most valence options first.
 */
function adjustForCharge(atoms: RawElement[], charge: number): KekuleElement[] {
  const adjusted = atoms.map((a) => ({ ...a, valence: getValence(a.name) }));
  if (charge === 0) return adjusted;

  let delta = -charge;

  // Try to absorb delta by switching atoms to alternate valence options
  for (const atom of adjusted) {
    if (delta === 0) break;
    const options = getValenceOptions(atom.name).sort((a, b) =>
      Math.abs(a - atom.valence + delta) - Math.abs(b - atom.valence + delta)
    );
    // Pick the option closest to current + delta
    const target = atom.valence + delta;
    const best = options.reduce((prev, cur) =>
      Math.abs(cur - target) < Math.abs(prev - target) ? cur : prev
    );
    delta -= best - atom.valence;
    atom.valence = best;
  }

  return adjusted;
}

/** Check the sum of all valences is even. */
function totalValenceIsEven(atoms: KekuleElement[]): boolean {
  const sum = atoms.reduce((acc, a) => acc + a.valence, 0);
  return sum % 2 === 0;
}

function tryCalculateKekuleStructure(molecule: Molecule, elements: RawElement[]): LewisStructure | undefined {
   // 2. Adjust valences for overall charge
  const atoms = adjustForCharge(elements, molecule.charge);

  // 3. All valences must be non-negative
  if (atoms.some((a) => a.valence < 0))
    return undefined;

  console.log(`${molecule.text}: Passed negative valence`);
  // 4. Total valence must be even
  if (!totalValenceIsEven(atoms))
    return undefined;

  console.log(`${molecule.text}: Even valence`);
  const kekule = getKekuleStructure(molecule);
  if (kekule === undefined)
    return undefined;

  console.log(`${molecule.text}: Valid kekule`);
  return undefined;
}

////////////////////////////////////////////////////////////////////////////////
// Lewis structure
////////////////////////////////////////////////////////////////////////////////

/**
 * Determine the central atom with the following rules:
 * 1. Single element: that atom is the center
 * 2. Carbon is always central in organic molecules
 * 3. Least electronegative non-hydrogen atom
 */
function pickCentralAtom(elements: RawElement[]): ElementName {
  let unique = elements.map((e) => e.name);
  unique = [...new Set(unique)];

  if (unique.length === 1)
    return unique[0];

  // Remove H as candidate (H is always terminal)
  const nonH = unique.filter((s) => s !== "H");
  if (nonH.length <= 1)
    // If H2 or some other atom
    return nonH.length === 0 ? unique[0] : nonH[0];

  // Pick least electronegative among non-H
  return nonH.reduce((best, sym) => {
    const bEN = periodicTable[best].electronegativity ?? 999;
    const sEN = periodicTable[sym].electronegativity ?? 999;
    return sEN < bEN ? sym : best;
  });
}

/** Cached results from `estimateMaxBonds` */
const cachedMaxBonds = new Map<string, number>();
/** For estimating happy bonds */
const MIN_BONDS: ElementRecord<number> = {
  C: 4, N: 3, O: 2,
  F: 1, Cl: 1, Br: 1, I: 1
};

/** See `estimateMaxBonds` for info. */
function estimateMaxBondsNoCache(element: ElementName): number {
  const { period, numberOfValence } = periodicTable[element];
  const canExpandOctet = period >= 3;

  if (numberOfValence === null)
    throw new Error(`${element} has unknown valence`);

  // 1. Noble gases
  // Full valence shell (8 e-, or 2 for He) means no unpaired electrons
  if (numberOfValence === 8) {
    // Period 3+ noble gases (Xe, Kr) can form bonds via d-orbital expansion
    return canExpandOctet ? 6 : 0;
  } else if (numberOfValence === 2 && period === 1) {
    return 0; // He
  }

  // 2. Period 1 (H only at this point)
  // Hydrogen obeys the duet rule — 1 bond fills its 1s shell
  if (period === 1) return 1;

  // 3. Period 2, strict octet, no d-orbital expansion
  // Ground-state unpaired electrons from Hund's rule on [2s^2][2p^n]:
  //   group 1 (Li): 1 unpaired  | 1 bond
  //   group 2 (Be): 0 unpaired  | can still form 2 bonds via sp hybridization
  //   group 13 (B): 1 unpaired  | 3 bonds via sp2 (empty p promotes)
  //   group 14 (C): 2 unpaired  | 4 bonds via sp3
  //   group 15 (N): 3 unpaired  | 3 bonds (+ lone pair)
  //   group 16 (O): 2 unpaired  | 2 bonds (+ 2 lone pairs)
  //   group 17 (F): 1 unpaired  | 1 bond  (+ 3 lone pairs)
  if (period === 2) {
    // Electrons needed to reach full octet = bonding capacity
    // Each bond contributes 1 electron from this atom toward filling to 8
    //
    // No atom can bond more than 4 times in period 2
    // Note: Be (VE=2) -> min(4,6)=4 but practically forms 2; this is a known simplification.
    // Lewis theory doesn't capture sp hybridization promotion, 4 is the safe upper bound.
    return Math.min(4, 8 - numberOfValence);
  }

  // 4. Period 3+ with d-orbital expansion
  // Ground-state unpaired electrons (Hund's rule on ns²np^x + available nd):
  //   ve <= 4  -> valenceElectrons itself (e.g. Si ve=4 -> 4 bonds)
  //   ve == 5  -> 5 unpaired via d promotion (e.g. P: 3 ground -> 5 promoted)
  //   ve == 6  -> 6 unpaired via d promotion (e.g. S: 2 ground -> 4 or 6)
  //   ve == 7  -> 7 unpaired is not chemically achieved; max observed = 7
  //               but Cl practically stays at 1 normally, up to 7 theoretically
  if (canExpandOctet) {
    if (numberOfValence <= 4) return numberOfValence;
    // For ve 5–7: maximum bonds = valenceElectrons (all electrons unpairable
    // by promoting into empty d-orbitals), capped at 7 (no element bonds 8x)
    return Math.min(numberOfValence, 7);
  }

  // Fallback: period 3+ without d-orbital expansion (metals like Na, Mg)
  // These don't covalently expand, treat them like period 2
  return Math.min(4, 8 - numberOfValence);
}

/**
 * Estimates the maximum number of bonds an element can form,
 * derived from its electronic structure rather than hardcoded.
 *
 * Rules applied (in order):
 *  1. Noble gases (full valence shell, no unpaired e-) -> 0
 *     Exception: period 3+ noble gases like Xe can expand octet -> up to 6
 *  2. Period 1 (H, He) -> capped at 1 (duet rule)
 *  3. Period 2 non-metals -> capped at 4 (strict octet, no d-orbital expansion)
 *  4. Period 3+ -> can expand octet using d-orbitals:
 *     maxBonds = min(8 - valenceElectrons, available d-orbital slots)
 *     but at minimum equals the number of unpaired electrons in ground state
 */
function estimateMaxBonds(element: ElementName): number {
  let out = cachedMaxBonds.get(element); 
  if (out !== undefined) return out;
  out = estimateMaxBondsNoCache(element);
  cachedMaxBonds.set(element, out);
  return out;
}

function estimateHappyBonds(element: ElementName): number {
  if (element in MIN_BONDS)
    return MIN_BONDS[element]!;
  // Fallback
  return estimateMaxBonds(element);
}

function collectElementCounts(elements: RawElement[]): ElementRecord<number> {
  let out: ElementRecord<number> = {}
  for (const el of elements) {
    if (out[el.name])
      out[el.name]! += el.count;
    else
      out[el.name] = el.count;
  }
  return out;
}

// Try rings if central molecules >5-6

// L-5-hydroxytryptophan > C11H12N2O3

// Benzaldehyde > C7H6O / C6H5CHO > (aromatic)
// molybdenum hexacarbonyl > C6MoO6 / Mo(CO)6
// Hydrogen peroxide > H2O2 / H-O-O-H
// Acetophenone > C8H8O / C6H5COCH3
// Ferrous oxide > Fe2O3 / O(FeO)2 / O=Fe−O−Fe=O
// Perchloric acid > HClO4 / H(ClO4)
// Ketene > C2H2O / CH2CO / H2C=C=O
// Disulfur decafluoride > S2F10 / (SF5)2
// Dithionic acid > H2S2O6 / (SO2OH)2
// Mercury(I) chloride > Hg2Cl2 / Cl-Hg-Hg-Cl
// Diborane > B2H6 / H2=B=H2=B=H2

// Periodyl periodate > I2O7 / IOOOOOOOI

// Thiobenzoic acid > C7H6OS / C6H5COSH
//   C6H5C(O)Cl + KSH → C6H5C(O)SH + KCl

// Methanethiol > CH4S / CH3SH
// Dimethyl disulfide > C2H6S2 / CH3SSCH3 / S2(CH3)2
//   2CH3SH + I2 → CH3SSCH3 + 2HI

// Ammonium > NH4{+}

interface Atom2 {
  id: number;
  name: ElementName;
  bonds: Set<number>;
  bondsLeft: number;
};

class Atom3 {
  id: number;
  name: ElementName;
  readonly maxBonds: number;
  bonds: [] = [];

  constructor(id: number, name: ElementName) {
    this.id = id;
    this.name = name;
    this.maxBonds = estimateMaxBonds(name);
  }

  bondsLeft() {
    return this.maxBonds - this.bonds.length;
  }
};

function calculateMultiCenter(data: MoleculeData): LewisStructure | undefined {
  const centerType = data.central;
  const prefCenterBonds = estimateHappyBonds(centerType);

  const counts = new Map<ElementName, number>();
  for (const [k, v] of counts) {
    if (k !== centerType)
      counts.set(k, v);
  }

  //////////////////////////////////////////////////////////////////////////////
  // 1. Set up central atoms

  // Create central atoms
  let idx = data.counts[centerType]!;
  const atoms: Atom2[] = [];
  for (let i = 0; i < idx; ++i) {
    atoms.push({
      id: i, name: centerType,
      bonds: new Set<number>(),
      bondsLeft: prefCenterBonds
    });
  }

  function assertID(id: number) {
    if (id >= atoms.length)
      throw new Error(`ID(a) out of bounds: ${id}`);
  }

  function addBond(idA: number, idB: number): boolean {
    if (idA === idB)
      return false;
    assertID(idA); assertID(idB);
    const a = atoms[idA], b = atoms[idB];
    // Try bonding
    if (a.bonds.has(idB)) {
      if (!b.bonds.has(idA))
        throw new Error(`Invalid state: ${idA} bonded to ${idB} but not inverse`);
      return true;
    }
    // Check
    if (a.bondsLeft === 0 || b.bondsLeft === 0) {
      return false;
    }
    a.bonds.add(idB);
    --a.bondsLeft;
    b.bonds.add(idA);
    --b.bondsLeft;
    return true;
  }
  
  // Bond central atoms with each other?
  if (prefCenterBonds > 1) {
    // Bond all atoms with the one next to them.
    for (let i = 0; i < idx - 1; ++i)
      addBond(i, i + 1);
  }

  return undefined;
}

function calculateSingleCenter(data: MoleculeData): LewisStructure | undefined {
  const isMetal = periodicTable[data.central].metal;
  return undefined;
}

export function calculateLewisStructure3(molecule: Molecule): LewisStructure | undefined {
  let data = new MoleculeData(molecule);
  data.log(`centralSymbol ${data.central} (multi: ${data.multiCenter})`);
  if (data.multiCenter)
    return calculateMultiCenter(data);
  else
    return calculateSingleCenter(data);
}
