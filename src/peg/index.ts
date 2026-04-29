import { parse as _parse } from './parser.js';
export { SyntaxError } from './parser.js';
import { periodicTable } from '../elements';
import { WeakMappingCache, getUniqueID } from '../object-utils';

/** Valid element names */
export type ElementName = keyof typeof periodicTable;

/** Single element */
export interface RawElement {
  count: number,
  name: ElementName
}

/** Single element, tagged */
export interface Element {
  element: true,
  count: number,
  name: ElementName
};

/** Group of elements, tagged */
export interface MoleculePart {
  element: false,
  count: number,
  data: Element[]
};

/** Set of elements and groups */
export interface Molecule {
  coefficient: number,
  charge: number,
  hydrate: number,
  data: (Element | MoleculePart)[],
  text: string
};

/** One set of Molecules for each side */
export interface Redox {
  lhs: Molecule[],
  rhs: Molecule[]
};

/** Just the base molecule */
type SimpleMolecule = Omit<Molecule, 'coefficient' | 'text'>;
/** Flattened molecule */
export type FlatMolecule = Record<string, number>;

/** Generates a unique ID for each molecule type */
function hashMoleculeNoCache<T extends SimpleMolecule>(mol: T): [string, FlatMolecule] {
  // Flatten molecule
  const data: Record<string, number> = {};
  function addToData(k: string, v: number) {
    if (k in data)
      data[k] += v;
    else
      data[k] = v;
  };

  for (const elOrData of mol.data) {
    const count = elOrData.count;
    if (elOrData.element) {
      addToData(elOrData.name, count);
    } else for (const el of elOrData.data) {
      addToData(el.name, el.count * count);
    }
  }

  // Get a unique identifier
  const ID = getUniqueID(data);
  const { charge, hydrate } = mol;
  return [`${ID}{${charge}*${hydrate}`, data];
}

/** Stores weak references to previously hashed molecules */
const prevHashedMoleculeMap = new WeakMappingCache<SimpleMolecule, [string, FlatMolecule]>(hashMoleculeNoCache);

/** Gets the number of elements per molecule, caches per object */
export function countElements<T extends SimpleMolecule>(mol?: T): FlatMolecule {
  if (mol === undefined || mol.data.length === 0)
    return {};
  // Try cache
  return prevHashedMoleculeMap.lookup(mol)[1];
}

/** Generates a unique ID for each molecule type, caches per object */
export function hashMolecule<T extends SimpleMolecule>(mol?: T): string {
  if (mol === undefined || mol.data.length === 0)
    return '';
  // Try cache
  return prevHashedMoleculeMap.lookup(mol)[0];
}

export function parse(input: string): Redox {
  return _parse(input) as Redox;
}

export function parseMolecule(input: string): Molecule {
  return _parse(input, { startRule: "Molecule" }) as Molecule;
}
