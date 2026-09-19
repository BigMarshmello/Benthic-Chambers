"""Stage 0 tests: the parameter file is structurally sound and honest."""

import math

import pytest

from benthic import config


@pytest.fixture(scope="module")
def cfg():
    return config.load()


def test_config_parses_and_validates(cfg):
    problems = config.validate(cfg)
    assert problems == [], "config.yaml validation problems:\n" + "\n".join(problems)


def test_every_parameter_has_a_status(cfg):
    params = list(config.walk(cfg))
    assert len(params) > 50
    for dotted, node in params:
        assert node["status"] in config.VALID_STATUS, dotted


def test_numeric_parameters_declare_units(cfg):
    for dotted, node in config.walk(cfg):
        value = node.get("value")
        if isinstance(value, (int, float)) and not isinstance(value, bool):
            assert "units" in node, f"{dotted} has a numeric value but no units"


def test_placeholders_and_assumptions_cite_a_source(cfg):
    for dotted, node in config.walk(cfg):
        if node["status"] in {"PLACEHOLDER", "ASSUMPTION", "LITERATURE"}:
            assert node.get("source"), f"{dotted} has no source"


def test_get_refuses_unset_placeholder(cfg):
    # permeability has not been measured yet -- reading it must fail loudly
    # rather than silently handing back a made-up number.
    with pytest.raises(ValueError, match="unset PLACEHOLDER"):
        config.get(cfg, "sediment.permeability")


def test_unmeasured_parameters_are_the_expected_ones(cfg):
    unset = set(config.unset_placeholders(cfg))
    # These are the measurements the project still owes. If this list
    # changes, the stage 0 note must change with it.
    assert "sediment.permeability" in unset
    assert "sediment.median_grain_size" in unset
    assert "production.surface_irradiance" in unset


def test_chamber_geometry_is_flagged_as_uncertain(cfg):
    # The user could not remember the exact diameter; it must not be
    # presented as measured.
    assert cfg["chamber"]["diameter"]["status"] == "PLACEHOLDER"
    assert cfg["chamber"]["height_above_sediment"]["status"] == "PLACEHOLDER"


def test_derived_chamber_volume_is_not_hardcoded(cfg):
    assert cfg["chamber"]["water_volume"]["value"] is None
    assert cfg["chamber"]["water_volume"]["status"] == "DERIVED"


def test_light_variant_requires_production_term(cfg):
    # Guard the physics gap: a light chamber with production disabled is
    # not a model of a light chamber.
    assert cfg["chamber"]["variant"]["options"] == ["dark", "light"]
    assert "enabled" in cfg["production"]


def test_respiration_basis_is_explicit(cfg):
    # The phi-convention ambiguity that broke the brief's analytic formula
    # must be a declared switch, not an implicit choice.
    assert cfg["respiration"]["basis"]["value"] in {"porewater", "bulk"}


def test_stirrer_default_matches_huettel_gust_anchor(cfg):
    radius = cfg["chamber"]["diameter"]["value"] / 2.0
    dps = cfg["chamber"]["stirrer_pressure"]["value"]
    rim_gradient = 2.0 * dps / radius  # Pa/m
    assert math.isclose(rim_gradient, 20.0, rel_tol=1e-9)  # = 0.2 Pa/cm
