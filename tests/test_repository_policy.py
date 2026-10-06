"""Regression checks for bilingual agent guidance and portable documentation."""
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'tools'))
from check_repository import path_violations, permits_chinese, markdown_link_violations


class DocumentationPolicyTests(unittest.TestCase):
    def test_user_agent_guidance_can_be_bilingual(self):
        for name in ('AGENTS.md', 'src/AGENTS.md',
                     '.agents/skills/alod-terminology/SKILL.md',
                     '.agents/skills/alod-terminology/references/notation.md'):
            with self.subTest(name=name):
                self.assertTrue(permits_chinese(name))

    def test_ordinary_guides_and_evidence_keep_english_policy(self):
        for name in ('README.md', 'docs/guide.md', 'docs/provenance/data.json',
                     'ALOD_SUBPROJECT_AGENT_PLAN_20260909.md',
                     '.agents/skills/example/data.json',
                     'docs/.agents/skills/example/SKILL.md'):
            with self.subTest(name=name):
                self.assertFalse(permits_chinese(name))

    def test_portability_check_still_applies_to_bilingual_text(self):
        text = '\u8bf4\u660e: ' + '/' + 'home/example/project\nhttps://example.org/home/project\n'
        self.assertEqual(path_violations(text), [1])
        self.assertEqual(path_violations('\u8bf4\u660e: $PROJECT_ROOT/results'), [])

    def test_removed_guide_is_reported_but_external_links_are_not_local(self):
        text='[ok](runtime.md#recovery) [old](old.md) [web](https://example.org/old.md)'
        errors=markdown_link_violations('docs/README.md',text,lambda name: name=='docs/runtime.md')
        self.assertEqual(len(errors),1)
        self.assertIn('old.md',errors[0])

    def test_encoded_relative_links_are_resolved(self):
        errors=markdown_link_violations('docs/guide.md','[a](../data/a%20b.json)',
                                        lambda name: name=='data/a b.json')
        self.assertEqual(errors,[])

    def test_local_planning_file_cannot_supply_public_document_link(self):
        from check_repository import LOCAL_ONLY_PLAN
        errors=markdown_link_violations('README.md',f'[plan]({LOCAL_ONLY_PLAN})',
                                       lambda name: name in {'README.md','docs/runtime.md'})
        self.assertEqual(len(errors),1)


if __name__ == '__main__':
    unittest.main()
