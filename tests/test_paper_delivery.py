import sys
import unittest
import json
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from paper.data import adaptive, clean, native
from paper.selection import completed_cycles, target_prefix, promotions


class PaperDeliveryTests(unittest.TestCase):
    def fixture(self):
        return dict(
            config=dict(member_ids=[0, 1]),
            states=[dict(state_id=0, coarse_free=8, rank=2, eta=[0.3])],
            audits=[dict(state_id=0, sample=i, E=0.1, exact_norm=2, role="train") for i in (0, 1)],
        )

    def test_missing_and_duplicate_audits(self):
        x = self.fixture()
        x["audits"].pop()
        with self.assertRaises(ValueError):
            adaptive(x)
        x = self.fixture()
        x["audits"].append(x["audits"][0])
        with self.assertRaises(ValueError):
            adaptive(x)

    def test_no_training_audit_confusion(self):
        r = adaptive(self.fixture())
        self.assertNotIn("test_count", r["rows"][0])
        self.assertEqual(r["rows"][0]["N_on"], 10)

    def test_non_nominal_residual_gate(self):
        for residual in (1e-4, float("nan"), -1.0):
            data = self.fixture()
            data["audits"][1]["PG_residual"] = residual
            with self.assertRaises(ValueError):
                adaptive(data)

    def test_native_mesh_does_not_use_an_arbitrary_audit(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            (folder / "run.json").write_text(json.dumps(dict(config={}, status="completed")))
            (folder / "states.jsonl").write_text(
                json.dumps(dict(state=0, free_dof=8, relative_energy_error=0.1)) + "\n"
            )
            audit = folder / "audits/intermediate"
            audit.mkdir(parents=True)
            (audit / "mesh.json").write_text("{}")
            self.assertNotIn("mesh", native(folder))
            (folder / "mesh.json").write_text('{"terminal":true}')
            self.assertEqual(native(folder)["mesh"], {"terminal": True})

    def test_portable_config_whitelist(self):
        r = adaptive(self.fixture())
        r["config"]["host_workspace"] = "private-host-location"
        result = clean({"e1-ALOD": r})
        self.assertNotIn("host_workspace", result["e1-ALOD"]["config"])
        r["rows"][0]["E"] = float("nan")
        with self.assertRaises(ValueError):
            clean({"e1-ALOD": r})

    def test_e3_crossing_in_an_intermediate_sweep(self):
        run = dict(
            rows=[
                dict(state=0, N_on=8, E=0.1, cycle_complete=False),
                dict(state=1, N_on=8, E=0.009, cycle_complete=False),
                dict(state=2, N_on=12, E=0.008, cycle_complete=True),
            ]
        )
        self.assertEqual([r["state"] for r in target_prefix(run)], [0, 1])
        self.assertEqual([r["state"] for r in completed_cycles(run)], [0, 2])
        run["rows"][1]["E"] = 0.02
        run["rows"][2]["E"] = 0.015
        self.assertEqual(len(target_prefix(run)), 3)  # Preserve resource-stopped endpoints.

    def test_promotion_events_survive_portable_roundtrip(self):
        data = self.fixture()
        data["ell_checks"] = [
            dict(action="promote", state_id=0, ell=2, next_ell=3),
            dict(action="keep", state_id=0, ell=3, next_ell=3),
        ]
        run = clean({"e1-ALOD": adaptive(data)})["e1-ALOD"]
        self.assertEqual(promotions(run), [dict(state_id=0, ell=2, next_ell=3)])


if __name__ == "__main__":
    unittest.main()
