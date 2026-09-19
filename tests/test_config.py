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


def test_chamber_bore_is_now_measured(cfg):
    # Revision 1: the production chamber bore and OD are known.
    assert cfg["chamber"]["bore"]["status"] == "MEASURED"
    assert cfg["chamber"]["bore"]["value"] == pytest.approx(0.168)
    assert cfg["chamber"]["outer_diameter"]["value"] == pytest.approx(0.200)


def test_tube_length_is_still_uncertain(cfg):
    # "maybe around 60 cm, definitely less than 1 m" -- not a measurement.
    assert cfg["chamber"]["tube_length"]["status"] == "PLACEHOLDER"
    assert cfg["chamber"]["tube_length"]["range"] == [0.40, 1.00]


def test_geometry_derived_from_measurements_is_not_hardcoded(cfg):
    # Water column height and volume follow from tube length and
    # insertion depth; hardcoding them would let a stale number survive
    # a change to either input.
    for name in ("height_above_sediment", "water_volume"):
        assert cfg["chamber"][name]["value"] is None, name
        assert cfg["chamber"][name]["status"] == "DERIVED", name


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
    # The anchor is a pressure GRADIENT (0.2 Pa/cm), so delta_p_s must
    # track the bore. Revision 1 rescaled it from 1.50 Pa (R = 0.150 m)
    # to 0.84 Pa (R = 0.084 m).
    radius = cfg["chamber"]["bore"]["value"] / 2.0
    dps = cfg["chamber"]["stirrer_pressure"]["value"]
    rim_gradient = 2.0 * dps / radius  # Pa/m
    assert math.isclose(rim_gradient, 20.0, rel_tol=1e-2)  # = 0.2 Pa/cm


def test_stirrer_transfer_efficiency_is_flagged_as_uncalibrated(cfg):
    # The tall-tube swirl-decay problem. Default is the OPTIMISTIC bound,
    # so the model overstates rather than understates chamber advection
    # until it is calibrated.
    node = cfg["chamber"]["stirrer_transfer_efficiency"]
    assert node["status"] == "ASSUMPTION"
    assert node["value"] == 1.0
    assert node["range"][0] <= 0.2


def test_aspect_ratio_limit_reflects_published_chamber_designs(cfg):
    assert cfg["chamber"]["aspect_ratio_max"]["value"] == pytest.approx(2.0)
