with open("example.txt", "r") as infile:
    numbers = [line.split()[0] for line in infile]

numbers.sort(key=int)

with open("input.txt", "w") as outfile:
    outfile.write("\n".join(numbers) + "\n")
