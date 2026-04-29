{{
function makeNested(data, count) {
  return { element: false, count: count ?? 1, data: data };
}
}}

Redox "redox reaction"
  = lhs:Expression _ Arrow rhs:Expression {
  	return { lhs: lhs, rhs: rhs };
  }

Expression
  = head:Molecule _ "+" tail:Expression { return [head, ...tail]; }
  / data:Molecule { return [data]; }

Molecule "molecule"
  = coeff:Natural? _ data:MoleculePart+ hydrate:Hydrate? charge:Charge? {
    return {
      coefficient: coeff ?? 1,
      charge: charge ?? 0,
      hydrate: hydrate ?? 0,
      data: data,
      text: text().trim()
    };
  }

MoleculePart "part of top-level molecule"
  = _ "(" data:Element+ _ ")" count:Natural? { return makeNested(data, count); }
  / _ "[" data:Element+ _ "]" count:Natural? { return makeNested(data, count); }
  / Element

Element "element"
  = _ name:([A-Z][a-z]?) _ count:Natural? {
    return { element: true, count: count ?? 1, name: name.join('') };
  }

Hydrate
  = _ "*" value:Natural? _ "H" _ "2" _ "O" { return value ?? 1; }

Charge
  = _ "{" _ value:Integer? sign:("+" / "-")? _ "}" {
    value ??= 1;
    return sign !== "-" ? value : -value;
  }

Integer "integer"
  = _ (("-" / "+")? [0-9]+) { return parseInt(text(), 10); }

Natural "natural number"
  = _ [1-9][0-9]* { return parseInt(text(), 10); }

Arrow
  = "→" / "="+ / "-"* ">"

_ "whitespace"
  = [ \t\n\r]*
