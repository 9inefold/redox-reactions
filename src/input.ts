const MAX_TRIES = 3;

export async function getUserInput(question?: string): Promise<string> {
  var stdin = process.stdin, stdout = process.stdout;
  stdin.resume();
  if (question)
    stdout.write(question);

  return new Promise<string>((resolve, _) => {
    stdin.once('data', (buf: Buffer) => {
      resolve(buf.toString().trim());
    });
  });
}

export async function ask(question: string, format: RegExp, tries?: number) : Promise<string> {
  const maxTries = tries ?? MAX_TRIES;
  let failCount = 0;
  while (failCount < maxTries) {
    const data = await getUserInput(question);
    if (format.test(data))
      return data;
    // Invalid input, try again.
    ++failCount;
    process.stdout.write("Input should match: " + format + "\n");
  }
  
  throw new Error("Too many tries!");
}

type Transformer = (data: string) => string | undefined;

export async function askWith(question: string, validator: Transformer, tries?: number) : Promise<string> {
  const maxTries = tries ?? MAX_TRIES;
  let failCount = 0;
  while (failCount < maxTries) {
    const data = await getUserInput(question);
    const newData = validator(data);
    if (newData !== undefined)
      return data;
    // Invalid input, try again.
    ++failCount;
  }
  
  throw new Error("Too many tries!");
}
