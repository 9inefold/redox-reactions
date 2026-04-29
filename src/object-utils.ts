export type LaxRecord<K extends keyof any, V> = Partial<Record<K, V>>;

////////////////////////////////////////////////////////////////////////////////
// UniqueID

export type NameToNumber = { [k: string]: number; };

export function getUniqueID<KV extends NameToNumber>(data: KV): string {
  return Object.keys(data).sort().map(k => `${k}${data[k]!}`).join('');
}

////////////////////////////////////////////////////////////////////////////////
// WeakMappingCache

type Mapping<K, V> = (arg: K) => V;

/** Caches data based on object references */
export class WeakMappingCache<K extends WeakKey, V extends Exclude<unknown, undefined>> {
  /** Caches references */
  private cache: WeakMap<K, V>;
  private mapping: Mapping<K, V>;

  constructor(mapping: Mapping<K, V>) {
    this.cache = new WeakMap<K, V>();
    this.mapping = mapping;
  }

  lookup(k: K): V {
    // Get cached value, if possible
    const cached = this.cache.get(k);
    if (cached !== undefined)
      return cached;
    // Compute the value
    const unique = this.mapping(k);
    this.cache.set(k, unique);
    return unique;
  }

  update(k: K, v?: V) {
    if (v !== undefined) {
      this.cache.set(k, v);
      return;
    }
    // Compute the value
    const unique = this.mapping(k);
    this.cache.set(k, unique);
  }
};
