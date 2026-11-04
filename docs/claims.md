# Claims and evidence

The [claim register](claim-register.md) lists C01–C102, O01–O04 and P01–P12
in numerical order. The [paper-order index](claim-reading-order.md) uses the
same records. Stable claim IDs differ from printed theorem numbers.

This is a derived reader edition. Private operational history, authority
excerpts and obsolete workflow labels are omitted. Its editorial qualifications
preserve mathematical premises, refutations and open cases. It neither replaces
the retained original records nor constitutes a new proof review.

| Disposition | Meaning |
| --- | --- |
| RECORDED_PROOF_WITH_PREMISES | Recorded proof or certificate, subject to the stated premises |
| ADOPTED_SCOPED | Internally accepted within the named scope; no external endorsement implied |
| CONDITIONAL_EXECUTION | Historical execution or analytic-library/model premises remain explicit |
| REFUTED_AS_STATED | False as stated; follow the linked correction or replacement |
| SUPERSEDED_PARTIAL_RESULT | Replaced by a later result, with the later result's qualifications |
| SUPERSEDED_TRACKING | Earlier proposal or discovery retained for reference |
| HISTORICAL_PREPARATION | Preparatory work rather than a standalone current theorem |

The mathematical statements and proofs are in the papers. Internal argument
checks included model-assisted review; they are not external human specialist
endorsement. Hashes identify bytes, not mathematical truth or historical events.

## Important boundaries

- C18 is refuted; C31/C33 give the replacement common lattice and family.
- C69's earlier bound is sharpened by C70.
- C73's former B/C open cases are resolved only within C79's conditional
  classification. C78 records 107 opposite pairs per model.
- C89 leaves all 156 smaller-order Jacobian cases open.
- C49 is a retained generator extension, not an implicit premise of the shorter proof.
- C95 follows from the current C96 orbit argument.
- C102 retains lower bounds for the residual 20 classes, not exact minima.

Each direct reference has an exact digest and an explicit availability status
in [the register data](../data/claims/register.json). Omitted references are not
counted as verified. The [evidence index](../evidence/README.md) is a selection,
not a complete proof-dependency graph or runtime closure.

Run `python3 scripts/check-public-package.py` from the snapshot root for bounded
integrity and navigation checks. See [limitations](limitations.md) before
interpreting a successful check.
