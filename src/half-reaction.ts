import type { Redox, Molecule, ElementName } from './peg';
import { getOxidationStates, OxidationStates } from './oxidation';
import { formatCharge, formatOxidation, formatMolecule, formatReaction } from '.';

export interface OxidationDiff {
  reactant: Molecule,
  product: Molecule,
  element: ElementName,
  change: {
    from: number,
    to: number,
  }
};

export interface HalfReactions {
  oxidation: OxidationDiff,
  reduction: OxidationDiff,
};

function splitSymmetricDifference<T>(lhs: Set<T>, rhs: Set<T>): [T[], T[]] {
  return [[...lhs.difference(rhs)], [...rhs.difference(lhs)]];
}

export function getOxidationChange(diff: OxidationDiff): number {
  const { from, to } = diff.change;
  return to - from;
}

class PerElementRecord {
  /** Molecules containing this element */
  reverse: [Molecule, number][] = [];
  /** Set of seen states */
  seen = new Set<number>();

  add(mol: Molecule, val: number) {
    this.reverse.push([mol, val]);
    this.seen.add(val);
  }

  findFromOx(val: number): Molecule {
    const found = this.reverse.find(([_, v]) => val === v);
    if (found === undefined)
      throw new Error(`Could not find oxidation state ${formatCharge(val)}`);
    return found[0];
  }
};

class OxidationInfo {
  private molecules = new Map<Molecule, OxidationStates>();
  private elements = new Map<ElementName, PerElementRecord>();
  //private seen = new Set<ElementName>();

  add(mol: Molecule, ox: OxidationStates) {
    this.molecules.set(mol, ox);
    // Add reverse mappings
    for (const [el, val] of Object.entries(ox) as [ElementName, number][]) {
      let record = this.elements.get(el);
      // Init record if doesn't exist.
      if (record === undefined) {
        record = new PerElementRecord();
        this.elements.set(el, record);
      }
      // Register
      record.add(mol, val);
    }
  }

  findOxidationChanges(other: OxidationInfo): OxidationDiff[] {
    const out: OxidationDiff[] = [];
    const common = this.getCommonElts(other);
    // Search for differences
    for (const el of common) {
      const thisEl = this.elements.get(el)!,
            otherEl = other.elements.get(el)!;
      // Get the differences
      const [rDiff, pDiff] = splitSymmetricDifference(thisEl.seen, otherEl.seen);
      const thisLen = rDiff.length;
      // Check there aren't any complex setups
      if (thisLen !== pDiff.length)
        throw new Error(`Mismatched oxidation changes: ${thisLen} != ${pDiff.length}`);
      // No changes
      if (thisLen === 0)
        continue;
      // Ignore these for now
      if (thisLen > 1)
        throw new Error(`Multiple oxidation changes: ${thisLen} > 1`);
      // Get the reverse mappings
      const reactant = thisEl.findFromOx(rDiff[0]);
      const product = otherEl.findFromOx(pDiff[0]);
      // Add new entry
      out.push({
        reactant: reactant,
        product: product,
        element: el,
        change: {
          from: rDiff[0],
          to: pDiff[0],
        }
      });
    }
    return out;
  }

  private static readonly OH_SET = new Set<ElementName>(['O', 'H']);

  getCommonElts(other: OxidationInfo): ElementName[] {
    const thisKeys = this.getSeenElts(),
          otherKeys = other.getSeenElts();
    // Check only differences are in O/H
    const different = thisKeys
      .symmetricDifference(otherKeys)
      .difference(OxidationInfo.OH_SET);
    if (different.size !== 0) {
      const diffElts = [...different].join(', ');
      throw new Error(`Elements not balanced: ${diffElts}`);
    }
    // Get the intersecting elements
    const intersecting = thisKeys.intersection(otherKeys);
    return [...intersecting];
  }

  private getSeenElts(): Set<ElementName> {
    return new Set<ElementName>(this.elements.keys());
  }
};

/** Gets the oxidation states of each molecule on each side. */
function mapOxidationStates(equation: Redox): [OxidationInfo, OxidationInfo] {
  const unsupported: Molecule[] = [];

  function handleMolecules(mols: Molecule[]): OxidationInfo {
    const info = new OxidationInfo();
    // Get oxidation states for all molecules...
    for (const mol of mols) {
      const ox = getOxidationStates(mol);
      console.log(`${formatMolecule(mol)}: ${formatOxidation(ox)}`);
      // Unsupported oxidation state...
      if (ox === undefined) {
        unsupported.push(mol);
        continue;
      }
      info.add(mol, ox);
    }
    return info;
  }

  const lhsMols = handleMolecules(equation.lhs);
  console.log(' → ');
  const rhsMols = handleMolecules(equation.rhs);

  if (unsupported.length > 0) {
    const ls = unsupported.map(m => formatMolecule(m)).join(', ');
    throw new Error(`Sorry, these molecules are unsupported: ${ls}`);
  }

  console.log();
  return [lhsMols, rhsMols];
}

export function determineHalfReactions(equation: Redox): HalfReactions | undefined {
  const [lhs, rhs] = mapOxidationStates(equation);
  const changes: OxidationDiff[] = lhs.findOxidationChanges(rhs);
  const nchanges = changes.length;
  // Check for the right amount of changes.
  if (nchanges > 2)
    throw new Error(`Too many changes: ${nchanges} > expected (2)`);
  else if (nchanges < 2) {
    // Is this even valid? Sanity check.
    if (nchanges === 1)
      throw new Error('What? 1 change?');
    // Lel!
    console.log(`No changes needed for: ${formatReaction(equation)}`);
    return;
  }

  // Now we can find our reactions:
  // Oxidation: losing electrons (ox up)
  if (getOxidationChange(changes[0]) > 0)
    return { oxidation: changes[0], reduction: changes[1] };
  // Reduction: gaining electrons (ox down)
  else
    return { oxidation: changes[1], reduction: changes[0] };
}
