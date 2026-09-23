"""Exact product-row cancellation rates for the reduced Ring-LPN models."""

from __future__ import annotations

import csv
from pathlib import Path

from rev_cuckoo_full_lpn_field import (
    make_model,
    product_active_list,
    product_raw_addresses,
)
from rev_cuckoo_ring_lpn_channel_sweep import CONFIGS


def main() -> None:
    rows = []
    for config in CONFIGS:
        model = make_model(2, config.weight, config.block_length, config.modulus)
        observations = 0
        raw_rows = 0
        active_rows = 0
        cancelled_lists = 0
        for known_poly in model.poly_states:
            for hidden_poly in model.poly_states:
                for output_block in range(model.weight):
                    raw = product_raw_addresses(
                        model, known_poly, hidden_poly, output_block
                    )
                    active = product_active_list(
                        model, known_poly, hidden_poly, output_block
                    )
                    observations += 1
                    raw_rows += len(raw)
                    active_rows += len(active)
                    cancelled_lists += len(active) < len(raw)
        row = {
            "name": config.name,
            "modulus": config.modulus,
            "weight": config.weight,
            "block_length": config.block_length,
            "observations": observations,
            "mean_rows_after_address_dedup": raw_rows / observations,
            "mean_active_rows_after_value_dedup": active_rows / observations,
            "list_cancellation_probability": cancelled_lists / observations,
        }
        rows.append(row)
        print(
            f"{config.name}: raw={row['mean_rows_after_address_dedup']:.6f} "
            f"active={row['mean_active_rows_after_value_dedup']:.6f} "
            f"Pr[cancel]={row['list_cancellation_probability']:.6f}"
        )

    output = Path(__file__).with_name(
        "rev_cuckoo_product_cancellation_results.csv"
    )
    with output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=tuple(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {output}")


if __name__ == "__main__":
    main()
