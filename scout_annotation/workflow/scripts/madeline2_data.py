from pathlib import Path
import csv


class Sex:
    MALE = "M"
    FEMALE = "F"
    UNKNOWN = "."

    @classmethod
    def from_int(cls, sex):
        if sex == 1:
            return cls.MALE
        elif sex == 2:
            return cls.FEMALE
        else:
            return cls.UNKNOWN


class Status:
    AFFECTED = "A"
    UNAFFECTED = "U"
    UNKNOWN = "."

    @classmethod
    def from_int(cls, status):
        if status == 1:
            return cls.UNAFFECTED
        elif status == 2:
            return cls.AFFECTED
        elif status == 0 or status == -9:
            return cls.UNKNOWN
        else:
            raise NotImplementedError("continous phenotype values not supported")


class PedIndividual:
    def __init__(self, sample, father, mother, sex, status):
        self.sample = sample
        self.mother = mother if mother != "0" else "."
        self.father = father if father != "0" else "."
        self.sex = Sex.from_int(sex)
        self.status = Status.from_int(status)

    def __str__(self):
        return "\t".join([self.sample, self.sex, self.father, self.mother, self.status])


class Family:
    def __init__(self, family, individuals=None):
        self.family = family
        self.individuals = []
        if individuals is not None:
            self.individuals = individuals

    def add(self, individual):
        self.individuals.append(individual)

    def __repr__(self):
        return f"<Family {self.family} with {len(self)} individuals>"

    def __str__(self):
        lines = []
        for individual in self.individuals:
            lines.append(f"{self.family}\t{individual}")
        return "\n".join(lines)

    def __len__(self):
        return len(self.individuals)


def parse_ped(ped_file):
    families = {}
    with open(ped_file) as f:
        reader = csv.reader(f, delimiter="\t")
        for line in reader:
            if line[0] not in families:
                families[line[0]] = Family(line[0])
            ind = PedIndividual(line[1], line[2], line[3], int(line[4]), int(line[5]))
            families[line[0]].add(ind)
    return families


def madeline2_data(families):
    lines = [
        "\t".join(
            ["FamilyID", "IndividualID", "Gender", "Father", "Mother", "Affected"]
        )
    ]
    for fam in families.values():
        lines.append(str(fam))
    return "\n".join(lines)


def main():
    ped_file = Path(snakemake.input.ped)
    out_file = Path(snakemake.output.madeline2_data)

    families = parse_ped(ped_file)
    tsv_data = madeline2_data(families)

    with open(out_file, "w") as f:
        f.write(tsv_data)


if __name__ == "__main__":
    main()
