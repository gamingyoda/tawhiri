#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
tawhiri_json_to_csv.py

Tawhiri API の予測 JSON をフラットな CSV に変換します。

想定入力（例）:
- top-level に "prediction": [{"stage": "...", "trajectory":[{datetime, latitude, longitude, altitude}, ...]}, ...]
- ほかに "request", "metadata", "warnings" が入ることがあります

出力CSV（1ファイル）:
stage, point_index, datetime_utc, unix_time, seconds_from_launch, latitude, longitude, altitude_m

使い方:
  python tawhiri_json_to_csv.py pred.json out.csv
  python tawhiri_json_to_csv.py pred.json out.csv --sort-by-time
"""

from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple


def _parse_iso8601(dt_str: str) -> datetime:
    """ISO8601 (末尾Z/小数秒あり) を UTC aware datetime に変換"""
    s = dt_str.strip()
    if s.endswith("Z"):
        s = s[:-1] + "+00:00"
    # datetime.fromisoformat は 2025-...+00:00 / 小数秒OK
    dt = datetime.fromisoformat(s)
    if dt.tzinfo is None:
        # tz-naive は UTC とみなす
        dt = dt.replace(tzinfo=timezone.utc)
    else:
        dt = dt.astimezone(timezone.utc)
    return dt


def flatten_tawhiri_prediction(data: Dict[str, Any]) -> Tuple[List[Dict[str, Any]], datetime]:
    """
    予測JSONをフラットな行にして返す。
    戻り値:
      (rows, launch_dt)
    """
    pred = data.get("prediction")
    if not isinstance(pred, list) or len(pred) == 0:
        raise ValueError("JSON に 'prediction' 配列が見つかりません")

    # launch_datetime: request があればそれを優先。なければ最初の点。
    launch_dt: Optional[datetime] = None
    req = data.get("request")
    if isinstance(req, dict) and isinstance(req.get("launch_datetime"), str):
        launch_dt = _parse_iso8601(req["launch_datetime"])

    if launch_dt is None:
        # prediction[0].trajectory[0].datetime
        for block in pred:
            traj = block.get("trajectory")
            if isinstance(traj, list) and traj:
                d0 = traj[0].get("datetime")
                if isinstance(d0, str):
                    launch_dt = _parse_iso8601(d0)
                    break
    if launch_dt is None:
        raise ValueError("launch_datetime を推定できませんでした（request も trajectory も見つからない）")

    rows: List[Dict[str, Any]] = []
    for block in pred:
        stage = block.get("stage", "")
        traj = block.get("trajectory", [])
        if not isinstance(traj, list):
            continue
        for i, p in enumerate(traj):
            if not isinstance(p, dict):
                continue
            dt_str = p.get("datetime")
            lat = p.get("latitude")
            lon = p.get("longitude")
            alt = p.get("altitude")
            if not isinstance(dt_str, str):
                continue

            dt = _parse_iso8601(dt_str)
            unix_time = dt.timestamp()
            sec_from_launch = (dt - launch_dt).total_seconds()

            rows.append({
                "stage": stage,
                "point_index": i,
                "datetime_utc": dt.isoformat().replace("+00:00", "Z"),
                "unix_time": unix_time,
                "seconds_from_launch": sec_from_launch,
                "latitude": lat,
                "longitude": lon,
                "altitude_m": alt,
            })

    if not rows:
        raise ValueError("trajectory から有効な点を抽出できませんでした")

    return rows, launch_dt


def write_csv(rows: List[Dict[str, Any]], out_path: Path) -> None:
    fieldnames = [
        "stage",
        "point_index",
        "datetime_utc",
        "unix_time",
        "seconds_from_launch",
        "latitude",
        "longitude",
        "altitude_m",
    ]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)


def main() -> None:
    ap = argparse.ArgumentParser(description="Convert Tawhiri prediction JSON to CSV")
    ap.add_argument("input_json", type=Path, help="Tawhiri の予測 JSON ファイル")
    ap.add_argument("output_csv", type=Path, help="出力 CSV ファイル")
    ap.add_argument("--sort-by-time", action="store_true", help="datetime で昇順にソートして出力")

    args = ap.parse_args()

    data = json.loads(args.input_json.read_text(encoding="utf-8"))
    rows, launch_dt = flatten_tawhiri_prediction(data)

    if args.sort_by_time:
        rows.sort(key=lambda r: r["unix_time"])

    write_csv(rows, args.output_csv)

    print(f"[OK] wrote: {args.output_csv}")
    print(f"     points: {len(rows)}")
    print(f"     launch : {launch_dt.isoformat().replace('+00:00','Z')}")


if __name__ == "__main__":
    main()
