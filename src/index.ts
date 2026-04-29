import { getUserInput } from './input';
import { SyntaxError, parse, parseMolecule, hashMolecule, countElements } from './peg';
import type { Redox, MoleculePart, Molecule, Element, ElementName } from './peg';
import { determineHalfReactions, getOxidationChange } from './half-reaction'
import type { OxidationDiff, HalfReactions } from './half-reaction'
import { gaussianElimination } from './gauss';
import { getOxidationStates, OxidationStates } from './oxidation';

const CHOICES = "Enter 'e' for an equation, 't' to run tests, '?' for help, or 'q' to quit: ";
const CHOICES_REGEX = /[teh?q]/;

const EQUATION = "Enter a redox reaction: \n";

const SUPERSCRIPT = "⁰¹²³⁴⁵⁶⁷⁸⁹".split('');
const SUBSCRIPT = "₀₁₂₃₄₅₆₇₈₉".split('');

const H2O = "H₂O";

const CONST_H: Molecule = {
  coefficient: 1,
  charge: +1,
  hydrate: 0,
  data: [{ element: true, count: 1, name: 'H' }],
  text: 'H{+}'
};
const CONST_OH: Molecule = {
  coefficient: 1,
  charge: -1,
  hydrate: 0,
  data: [
    { element: true, count: 1, name: 'H' },
    { element: true, count: 1, name: 'O' }
  ],
  text: 'OH{-}'
};
const CONST_H2O: Molecule = {
  coefficient: 1,
  charge: 0,
  hydrate: 0,
  data: [
    { element: true, count: 2, name: 'H' },
    { element: true, count: 1, name: 'O' }
  ],
  text: 'H2O'
};

function chargeToSuperscript(charge: number): string {
  charge = Math.round(charge);
  // Ensure there is actually a charge number
  if (charge === 0)
    return '';

  // n?[+/-]
  const postfix = charge > 0 ? '⁺' : '⁻';
  charge = Math.abs(charge);
  // Ignore the number
  if (charge === 1)
    return postfix;

  // Parse it out
  let out = '';
  while (charge > 0) {
    out = SUPERSCRIPT[charge % 10] + out;
    charge = Math.floor(charge / 10);
  }
  return out + postfix;
}

function countToSubscript(count: number): string {
  count = Math.round(count);
  // Ignore the number
  if (count <= 1) {
    // Warn if <0?
    return '';
  }

  // Parse it out
  let out = '';
  while (count > 0) {
    out = SUBSCRIPT[count % 10] + out;
    count = Math.floor(count / 10);
  }
  return out;
}

export function formatMolecule(mol: Molecule): string {
  let out = "";

  function formatElement(el: Element) {
    return el.name + countToSubscript(el.count);
  };
  function formatMolpart(part: MoleculePart) {
    const els = part.data.map(el => formatElement(el)).join('');
    return `(${els})${countToSubscript(part.count)}`;
  };

  for (const elOrPart of mol.data) {
    if (elOrPart.element)
      out += formatElement(elOrPart);
    else
      out += formatMolpart(elOrPart);
  }

  if (mol.hydrate > 0) {
    out += ' ⋅ ';
    if (mol.hydrate === 1)
      out += H2O;
    else
      out += `${mol.hydrate}${H2O}`;
  }

  if (mol.charge)
    out = `[${out}]${chargeToSuperscript(mol.charge)}`;

  if (mol.coefficient > 1)
    out = `${mol.coefficient}${out}`;

  return out;
}

export function formatReaction({ lhs, rhs }: Redox): string {
  const reactant = lhs.map(m => formatMolecule(m)).join(' + ');
  const product = rhs.map(m => formatMolecule(m)).join(' + ');
  return reactant + ' → ' + product;
}

/** Format an individual charge */
export function formatCharge(state: number): string {
  return state > 0 ? `+${state}` : state.toString(10);
}
export function formatOxidation(ox: OxidationStates | undefined): string {
  if (ox === undefined) return '?';
  // Map to name/value pairs, map to string
  const entries = Object.entries(ox).map(([name, state]) => {
    if (Array.isArray(state))
      return `${name}: ...`;
    // Single entry
    return `${name}: ${formatCharge(state)}`;
  }).join(', ');
  return `{ ${entries} }`;
}

function formatHalfReaction({ reactant, product, element, change }: OxidationDiff): string {
  const { from, to } = change;
  const freact = `${formatMolecule(reactant)} → ${formatMolecule(product)}`;
  return `${element}[${formatCharge(from)} → ${formatCharge(to)}] ${freact}`;
}
function printHalfReactions(hr: HalfReactions) {
  const ox = formatHalfReaction(hr.oxidation);
  const rd = formatHalfReaction(hr.reduction);
  console.log(`Oxidation:\n  ${ox}\nReduction:\n  ${rd}`);
}

// Mn[(O2)2]2{3-} + FeSO4*7H2O -> MnO8*7H2O{-3} + FeSO4

async function getEquation(): Promise<Redox> {
  while (true) {
    const equation = await getUserInput(EQUATION);
    // Check if asking for help...
    if (equation === '?') {
      help();
      continue;
    }
    // Parse, or handle parse errors!
    try {
      return parse(equation);
    } catch (e) {
      if (e instanceof SyntaxError)
        console.error('Error: ' + e.message + '\n');
      else
        console.error('Error: ', e, '\n');
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

// From https://gist.github.com/bellbind/5468385accdee9df0d88
function gcd(a: number, b: number): number {
  // fast GCD aka Binary GCD
  if (a === 0) return b;
  if (b === 0) return a;
  if (a === b) return a;
  // remove even divisors
  var sa = 0;
  while (!(a & 1)) sa++, a >>= 1;
  var sb = 0;
  while (!(b & 1)) sb++, b >>= 1;
  var p = sa < sb ? sa : sb; // Power part of 2^p Common Divisor
  // euclidean algorithm: limited only odd numbers
  while (a !== b) {// both a and b should be odd
    if (b > a) {var t = a; a = b; b = t;} // swap as a > b
    a -= b; // a is even because of odd - odd
    do a >>= 1; while (!(a & 1)); // a become odd
  }
  return a << p; // Odd-Common-Divisor * 2^p
}

function coeffBalanceFactor(a: number, b: number): [number, number] {
  const C = gcd(a, b);
  return [b / C, a / C];
}

////////////////////////////////////////////////////////////////////////////////

function getElementAmount(mol: Molecule, name: ElementName): number {
  const counts = countElements(mol);
  return ((counts[name]) ?? 0) * mol.coefficient;
}

function balanceHalfReactions({ oxidation, reduction }: HalfReactions) {
  // Returns the number of electrons
  function balanceHR(diff: OxidationDiff): number {
    const { reactant, product, element } = diff;
    const get = (m: Molecule) => getElementAmount(m, element);
    let rv = get(reactant), pv = get(product);
    if (rv !== pv) {
      const C = gcd(rv, pv);
      reactant.coefficient *= (pv / C);
      product.coefficient *= (rv / C);
    }
    // Calculate electrons
    const e = Math.abs(getOxidationChange(diff));
    return e * get(reactant);
  }

  function getMulFactor(): [number, number] {
    // Balance the molecules so redox atoms are equal
    const oxE = balanceHR(oxidation);
    const rdE = balanceHR(reduction);
    return coeffBalanceFactor(oxE, rdE);
  }

  function mulHR(diff: OxidationDiff, C: number) {
    if (C === 1) return;
    const { reactant, product } = diff;
    reactant.coefficient *= C;
    product.coefficient *= C;
  }

  const [oxC, rdC] = getMulFactor();
  mulHR(oxidation, oxC);
  mulHR(reduction, rdC);
}

function balanceGaussian({ lhs, rhs }: Redox) {
  //const N = lhs.length + rhs.length + 1;
  const N = lhs.length + rhs.length;
  const mapping = new Map<string, number>();
  const matrix: number[][] = [];
  let id = 0;

  function makeRow(name: string): number[] {
    mapping.set(name, id++);
    const row = Array<number>(N).fill(0);
    matrix.push(row);
    return row;
  }

  function getRow(name: string): number[] {
    const off = mapping.get(name);
    if (off !== undefined)
      return matrix[off];
    return makeRow(name);
  }

  let off = 0;
  function fillColumn(mol: Molecule, C: number) {
    const counts = countElements(mol);
    for (const [name, val] of Object.entries(counts)) {
      //if (name === 'H' || name === 'O') continue;
      if (name === 'H') continue;
      const row = getRow(name);
      row[off] = val * C;
    }
    ++off;
  }

  for (const mol of lhs) {
    const { coefficient } = mol;
    fillColumn(mol, coefficient);
  }
  // Same thing, but it has to be negative
  for (const mol of rhs) {
    const { coefficient } = mol;
    fillColumn(mol, -coefficient);
  }

  const joined = [...lhs, ...rhs];
  const unknowns = Array<number>(matrix.length).fill(0);
  //console.log(matrix);
  
  const results = gaussianElimination(matrix, unknowns).map(v => -v);
  if (results.filter(v => Number.isNaN(v)).length > 0) {
    console.log('Error with gaussian elimination.');
    return;
  }

  //console.log(matrix);
  //console.log(results);

  for (const [off, val] of results.map(v => Math.round(v)).entries()) {
    if (val < 1) continue;
    joined[off].coefficient *= val;
  }
}

function sumCharges(mols: Molecule[]): number {
  let sum = 0;
  for (const {charge, coefficient} of mols) {
    sum += charge * coefficient;
  }
  return sum;
}

function addAcidic(equation: Redox) {
  const L = sumCharges(equation.lhs);
  const R = sumCharges(equation.rhs);
  // Make a new H
  let H = { ...CONST_H };
  H.coefficient = Math.abs(L - R);
  if (L < R) equation.lhs.push(H);
  else if (L > R) equation.rhs.push(H);
}

function addBasic(equation: Redox) {
  const L = sumCharges(equation.lhs);
  const R = sumCharges(equation.rhs);
  // Make a new H
  let OH = { ...CONST_OH };
  OH.coefficient = Math.abs(L - R);
  if (L > R) equation.lhs.push(OH);
  else if (L < R) equation.rhs.push(OH);
}

function balanceOxygen(equation: Redox) {
  function sumOxygen(mols: Molecule[]): number {
    let sum = 0;
    for (const mol of mols)
      sum += getElementAmount(mol, 'O');
    return sum;
  };
  const L = sumOxygen(equation.lhs);
  const R = sumOxygen(equation.rhs);
  // Make a new H
  let H2O = { ...CONST_H2O };
  H2O.coefficient = Math.abs(L - R);
  if (L > R) equation.rhs.push(H2O);
  else if (L < R) equation.lhs.push(H2O);
}

function checkBalance({ lhs, rhs }: Redox): boolean {
  const m = new Map<string, [number, number]>();
  function get(name: string): [number, number] {
    if (!m.has(name))
      m.set(name, [0, 0]);
    return m.get(name)!;
  }
  function addUp(mol: Molecule[], idx: 0 | 1) {
    // Add it all up
    for (const m of mol) {
      const counts = countElements(m);
      for (const [k, v] of Object.entries(counts))
        get(k)[idx] += v * m.coefficient;
    }
  }
  // do it
  addUp(lhs, 0);
  addUp(rhs, 1);
  // Check
  let out = true;
  for (const [k, [l, r]] of m.entries()) {
    if (l !== r) {
      out = false;
      console.log(`${k}: ${l} != ${r}`);
    }
  }
  return out;
}

function handleEnd(equation: Redox, type?: string) {
  balanceOxygen(equation);
  //balanceGaussian(equation);
  const balanced = checkBalance(equation);

  console.log(`\n${type ?? 'Result'}:`);
  console.log(' ', formatReaction(equation));
  console.log(' ', balanced ? 'balanced!' : 'NOT balanced...');
}

function handleRedox(equation: Redox) {
  const HR = determineHalfReactions(equation);
  if (HR === undefined) return;
  printHalfReactions(HR);
  balanceHalfReactions(HR);
  //printHalfReactions(HR);
  balanceGaussian(equation);

  const equationCopy = {
    lhs: [...equation.lhs],
    rhs: [...equation.rhs]
  };

  addAcidic(equation);
  handleEnd(equation, 'Acidic');
  addBasic(equationCopy);
  handleEnd(equationCopy, 'Basic');
}

async function handleEquation() {
  //console.log('Equation');
  const equation = await getEquation();
  //console.log("\nParsing equation...");
  handleRedox(equation);
  console.log();
}

function runTests() {
  function test(eq: string) {
    console.log('\x1b[36m%s\x1b[0m', `${eq}: `);
    try {
      handleRedox(parse(eq));
    } catch (e) {
      console.error(e);
    }
    console.log();
  };

  // Cr2O7^2- + HNO2 -> Cr^3+ + NO3^-

  test('Cr2O7{2-} + HNO2 -> Cr{3+} + NO3{-}');
  test('Cu + NO3{-} -> Cu{2+} + NO2');
  test('Cu + 2Ag{+} -> Cu{2+} + 2Ag');
  test('Cu + Ag{+} -> Cu{2+} + Ag');
  test('Cu{+} + Fe -> Fe{3+} + Cu');
  test('TiO2 + Ti -> TiO');
  test('Br2 + F{-} -> F2 + Br{-}');
  test('Mn + H2(SO4) -> Mn(SO4) + S + H2O');
  test('Al + H2O = Al2O3 + H');
  test('Ag + H2S + O2 -> Ag2S + H2O');
  
  //test('4Ag + 2H2S + O2 -> 2Ag2S + 2H2O');
  test('H2(SO4) -> 2H + SO4');
  quit(0);
}

async function main() {
  //printKnownLewis();

  function testMolecules(eq: string) {
    console.log(eq, ":");
    try {
      //handleEquation2(parse(eq));
    } catch (e) {
      console.error(e);
    }
    console.log();
  };

  //testMolecules('C2H5CHO + H2 → CH3CH2CH2OH');
  //testMolecules('2FeO(OH) -> 2FeO + H2O2');
  //testMolecules('CH3CH3 -> H2 + C2H4');
  //testMolecules('CH3COO{-} + Na{+} -> CH3COONa');
  //testMolecules('C6H5C(O)Cl + KSH -> C6H5C(O)SH + KCl')
  //testMolecules('CH3SH + I2 -> CH3SSCH3 + HI')

  //runTests();

  console.log('Welcome to the redox equation balancer!');
  while (true) {
    const choice = await getUserInput(CHOICES);
    if (!CHOICES_REGEX.test(choice)) {
      console.log(`Unknown command '${choice}', type '?' for help.`);
      continue;
    }

    console.error();
    switch (choice) {
    case 't':
      runTests();
      break;
    case 'e':
      await handleEquation();
      break;
    case '?':
    case 'h':
      help();
      break;
    case 'q':
      process.exitCode = 0;
      console.log('Quitting...');
      return;
    }
  }
}

function help() {
  console.log(`Syntax Help:

MOLECULES are written in the form:
  [COEFF] X[n]Yz[m] * [#HYDRATES]H2O {[CHARGE][+/-]}
For example:
  CH3COO{-}     -> Acetate
  Mg(ClO3)2     -> Magnesium Chlorate
  2FeSO4*7H2O   -> 2 x Iron(II) Sulfate Heptahydrate

Parentheses can be used to group ions, but only one level is currently supported.
For example, this is valid:
  (CH3)2
But this is not:
  (C(H)3)2

REACTIONS must be in the following form:
  nX + nY -> nZ
The arrow Regex is /→|=+|-*>/, meaning these are all valid:
  2CN{-} + Au ---> Au(CN)2{-}
  NH3 + ClO{-} > N2H4 + Cl{-}
  AlH4{-} + H2CO === Al{3+} + CH3OH
  Hg + ZnSO4 → HgO4S + Zn
`);
}

function quit(code?: number): void {
  process.exit(code);
}

////////////////////////////////////////////////////////////////////////////////

const KNOWN_LEWIS: [string, string][] = [
  ['OH{-}',     '{ H: +1, O: -2 }'], // Hydroxide
  ['H2O',       '{ H: +1, O: -2 }'], // Water
  ['H2O2',      '{ H: +1, O: -1 }'], // Peroxide
  // Sulf-
  ['SO4{2-}',   '{ O: -2, S: +6 }'], // Sulfate
  ['SO3{2-}',   '{ O: -2, S: +4 }'], // Sulfite
  // Carbon-
  ['CO3{2-}',   '{ C: +4, O: -2 }'], // Carbonate
  // Nitr-
  ['NO3{-}',    '{ N: +5, O: -2 }'], // Nitrate
  ['NO2{-}',    '{ N: +3, O: -2 }'], // Nitrite
  // Chlor-
  ['ClO4{-}',   '{ Cl: +7, O: -2 }'], // Perchlorate
  ['ClO3{-}',   '{ Cl: +5, O: -2 }'], // Chlorate
  ['ClO2{-}',   '{ Cl: +3, O: -2 }'], // Chlorite
  ['ClO{-}',    '{ Cl: +1, O: -2 }'], // Hypochlorite
  // Misc.
  ['CH3COO{-}', '{ C: +3, H: +1, O: -2 }'], // Acetate
];

function printKnownLewis() {
  const cache = new Set<string>();
  let mappings: [string, string][] = [];
  const invalids: string[] = [];
  const dupes: string[] = [];

  function handle([name, ox]: [string, string]) {
    try {
      const mol = parseMolecule(name);
      const hash = hashMolecule(mol);
      if (cache.has(hash)) {
        dupes.push(`${name}: ${ox}`);
        return;
      }
      mappings.push([hash, ox]);
      cache.add(hash);
    } catch (e) {
      invalids.push(`${name}: ${ox}`);
    }
  }

  for (const known of KNOWN_LEWIS)
    handle(known);

  const mapped = mappings.map(([name, ox]) => `'${name}': ${ox}`).join(',\n  ');
  console.log(`const KNOWN_STATES: Record<string, OxidationStates> = {\n  ${mapped}\n};\n`);

  if (dupes.length > 0)
    console.log("DUPE: ", dupes);
  if (invalids.length > 0)
    console.log("INVALID: ", invalids);
}

main().catch((reason) => {
  console.error(reason);
}).finally(() => {
  console.log('')
  quit(0);
});
