from __future__ import annotations

from .errors import JsonObject


def assign_profile_next_action() -> JsonObject:
    return {
        "action": "ask_user",
        "question": (
            "Give this genome a profile nickname (e.g. a first name or initials), and "
            "should it be the default profile for this machine?"
        ),
        "operation": "active_genome_index.assign_user_genome",
        "params": {
            "nickname": "<profile nickname>",
            "agi_id": "<active_genome_index.agi_id from this parse result>",
            "set_default_user": False,
        },
        "then": (
            "For an already-parsed genome, call active_genome_index.assign_user_genome "
            "with nickname and the parsed agi_id. During a new parse, pass "
            "user_nickname to genomi.parse_source instead."
        ),
    }


def reference_pass_next_action(job_id: object, job_path: object) -> JsonObject:
    action: JsonObject = {
        "action": "background_job",
        "operation": "active_genome_index.build_reference_pass",
        "why": (
            "Variants are ready now. The reference-block tail (~96% of a gVCF) is "
            "being appended in the background; coverage / 'is this site confirmed "
            "reference' answers stay provisional until it reports completed."
        ),
        "then": "Call genomi.check_background_job with this job_id to watch it finish.",
    }
    if job_id is not None:
        action["job_id"] = job_id
    if job_path is not None:
        action["job_path"] = job_path
    return action
