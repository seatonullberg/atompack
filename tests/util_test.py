from atompack.util import AttributeMap


def test_attribute_map():
    attmap = AttributeMap(name="test")
    assert attmap.name == "test"
    assert attmap["name"] == "test"
    for key, value in attmap.items():
        assert key == "name"
        assert value == "test"
    attmap["name"] = "newtest"
    assert attmap.name == "newtest"
    assert attmap["name"] == "newtest"
    assert len(attmap) == 1
    del attmap["name"]
    assert len(attmap) == 0
