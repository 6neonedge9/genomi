from __future__ import annotations

from typing import Any

from ...runtime.liftover import LiftoverConfigurationError, get_liftover

JsonObject = dict[str, Any]


def lift_score_variants(
    variants: list[JsonObject],
    *,
    source_build: str,
    target_build: str,
) -> dict[str, Any]:
    """Translate PGS-style score variants between genome builds.

    Lift each variant's (chrom, pos) using ``genomi.runtime.liftover``.
    Variants whose auto-generated ``variant_id`` was built from
    ``"<chrom>:<pos>:<effect>:<other>"`` are regenerated on the new
    coordinates so downstream logs/audits stay consistent. Records that
    fail to lift (chain gap, strand flip, missing coordinates) are
    returned in ``dropped`` with the reason recorded; this keeps the
    overall variant-accounting story honest without crashing the score.
    """

    try:
        lifter = get_liftover(source_build, target_build)
    except LiftoverConfigurationError:
        raise
    lifted: list[JsonObject] = []
    dropped: list[JsonObject] = []
    for variant in variants:
        chrom = variant.get("chrom")
        pos = variant.get("pos")
        if chrom is None or pos is None:
            dropped.append({**dict(variant), "liftover_reason": "missing_coordinates"})
            continue
        try:
            pos_int = int(pos)
        except (TypeError, ValueError):
            dropped.append({**dict(variant), "liftover_reason": "invalid_position"})
            continue
        result = lifter.lift_position_full(str(chrom), pos_int)
        if result is None:
            dropped.append({**dict(variant), "liftover_reason": "unmapped"})
            continue
        target_chrom, target_pos, strand = result
        if strand != "+":
            dropped.append({**dict(variant), "liftover_reason": "strand_flipped"})
            continue
        new_variant = dict(variant)
        new_variant["chrom"] = target_chrom
        new_variant["pos"] = target_pos
        new_variant["liftover"] = {
            "source_build": source_build,
            "target_build": target_build,
            "source_chrom": str(chrom),
            "source_pos": pos_int,
            "target_chrom": target_chrom,
            "target_pos": target_pos,
            "strand": strand,
        }
        old_auto_id = f"{chrom}:{pos_int}:{variant.get('effect_allele')}:{variant.get('other_allele') or ''}"
        if variant.get("variant_id") == old_auto_id:
            new_variant["variant_id"] = (
                f"{target_chrom}:{target_pos}:"
                f"{variant.get('effect_allele')}:{variant.get('other_allele') or ''}"
            )
        lifted.append(new_variant)
    return {
        "lifted": lifted,
        "dropped": dropped,
        "source_build": source_build,
        "target_build": target_build,
    }
