# Maintainer Checklist: Stale Branch Lifecycle

This document defines the **official maintainer procedure** for identifying, archiving,
and deleting stale branches in the UEDGE repository.

This checklist supports the documented policy that branches with no activity for an
extended period are preserved for history, made read-only, and eventually removed.

---

## Definition of a Stale Branch

A branch is considered **stale** if **all** of the following are true:

- No commits for **12 months**
- No open or merged pull requests for **12 months**
- Not referenced by:
  - an active release
  - an open issue or milestone
  - current development plans

---

## Phase 1: Identification

- [ ] Verify last commit date (≥ 12 months ago)
- [ ] Verify no PR activity in the last 12 months
- [ ] Confirm branch is not:
  - used for a release
  - referenced in active issues
  - required for ongoing work
- [ ] Identify likely original author(s) or owner(s)

---

## Phase 2: Notification

Before archiving, notify stakeholders.

- [ ] Open a tracking issue titled:
      **“Archiving stale branch: `<branch-name>`”**
- [ ] Include in the issue:
  - branch name
  - last commit date
  - reason for staleness
  - planned archive date
  - planned deletion date (archive + 12 months)
- [ ] Tag relevant maintainers and contributors
- [ ] Allow **2–4 weeks** for objections or retention requests

---

## Phase 3: Archive (Rename to `stale/*`)

Once the notification period has passed:

- [ ] Rename the branch to:
      ```
      stale/<original-branch-name>
      ```
- [ ] Verify branch protection is enabled:
  - no direct pushes
  - PRs blocked or restricted
  - force-push disabled
  - deletion disabled
- [ ] Add an entry to the stale branch log (see below)

> Archived branches are **read-only** and **not valid PR targets**.

---

## Phase 4: Documentation & Logging

Maintain a persistent record of archived branches.

- [ ] Add or update an entry in `docs/stale-branches-log.md` with:
  - original branch name
  - stale branch name
  - archive date
  - planned deletion date
  - link to tracking issue

Example entry:

```markdown
| Original Branch | Stale Branch | Archived | Delete After | Issue |
|-----------------|-------------|----------|--------------|-------|
| feature/old-x | stale/feature/old-x | 2024-06-01 | 2025-06-01 | #123 |
```

---

## Phase 5: Deletion (After Archive Period)

After **12 additional months** in the stale state:

- [ ] Reconfirm no renewed development need exists
- [ ] Post final notice comment on the tracking issue
- [ ] Delete the `stale/*` branch
- [ ] Update the stale branch log with the deletion date
- [ ] Close the tracking issue

---

## Phase 6: Revival (If Needed)

If development must resume:

- [ ] Create a new active branch from the stale branch tip:
      ```
      feature/<new-topic>
      feature/v9-<new-topic>
      ```
- [ ] Do **not** develop directly on `stale/*`
- [ ] Reference the original stale branch and issue in the new PR

---

## Notes & Best Practices

- Archiving is preferred over deletion to preserve scientific history
- Stale branches should **never** silently disappear
- This process ensures:
  - auditability
  - contributor respect
  - long-term repository hygiene

---

**This checklist is mandatory for all UEDGE maintainers.**
