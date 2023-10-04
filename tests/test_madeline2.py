from pathlib import Path
from pytest import fixture

from scout_annotation.workflow.scripts import madeline2_data


@fixture
def trio_ped():
    return Path(Path(__file__).parent, "integration", "data", "ceph1463_trio.ped")


@fixture
def empty_family():
    return madeline2_data.Family("fam01")


@fixture
def individuals():
    return lambda sex, status: madeline2_data.PedIndividual(
        "sample", "father", "mother", sex, status
    )


def test_family_id(empty_family):
    assert empty_family.family == "fam01"


def test_family_size(empty_family):
    assert len(empty_family) == 0


def test_individuals(individuals):
    ind = individuals(1, 1)
    assert ind.sex == madeline2_data.Sex.MALE == "M"
    assert ind.status == madeline2_data.Status.UNAFFECTED == "U"

    ind = individuals(2, 2)
    assert ind.sex == madeline2_data.Sex.FEMALE == "F"
    assert ind.status == madeline2_data.Status.AFFECTED == "A"

    ind = individuals(0, 0)
    assert ind.sex == madeline2_data.Sex.UNKNOWN == "."
    assert ind.status == madeline2_data.Status.UNKNOWN == "."


def test_parse_ped(trio_ped):
    families = madeline2_data.parse_ped(trio_ped)
    assert len(families) == 1
    assert len(families["ceph1463"]) == 3
    assert repr(families["ceph1463"]) == "<Family ceph1463 with 3 individuals>"


def test_madeline2_format(trio_ped):
    families = madeline2_data.parse_ped(trio_ped)
    tsv_data = madeline2_data.madeline2_data(families).splitlines()
    assert tsv_data[0].split("\t") == [
        "FamilyID",
        "IndividualID",
        "Gender",
        "Father",
        "Mother",
        "Affected",
    ]
    assert tsv_data[1].split("\t") == ["ceph1463", "NA12877", "M", ".", ".", "U"]
    assert tsv_data[2].split("\t") == ["ceph1463", "NA12878", "F", ".", ".", "U"]
    assert tsv_data[3].split("\t") == [
        "ceph1463",
        "NA12879",
        "F",
        "NA12877",
        "NA12878",
        "A",
    ]
