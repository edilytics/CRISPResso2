# Writers Module Refactor Design

## Overview

Move VCF output code from `CRISPRessoUtilities.py` to a new `writers/` module structure. This establishes a pattern for organizing output writers, with VCF as the first module.

## Directory Structure

**Before:**
```
CRISPResso2/
├── CRISPRessoUtilities.py
└── ...

tests/unit_tests/
└── test_CRISPRessoUtilities/
    ├── test_build_alt_map.py
    └── test_write_vcf.py
```

**After:**
```
CRISPResso2/
├── writers/
│   ├── __init__.py
│   └── vcf.py
└── ...

tests/unit_tests/
└── test_writers/
    └── test_vcf/
        ├── test_build_alt_map.py
        └── test_write_vcf.py
```

## Import Changes

**CRISPRessoCORE.py:**
```python
# Before
from CRISPResso2 import (CRISPRessoShared,
                         CRISPRessoPlot,
                         CRISPRessoUtilities)
# usage: CRISPRessoUtilities.write_vcf_file(...)

# After
from CRISPResso2 import (CRISPRessoShared,
                         CRISPRessoPlot)
from CRISPResso2.writers import vcf
# usage: vcf.write_vcf_file(...)
```

**Test files:**
```python
# Before
from CRISPResso2 import CRISPRessoUtilities as utilities
# usage: utilities.build_alt_map(...)

# After
from CRISPResso2.writers import vcf
# usage: vcf.build_alt_map(...)
```

## File Changes

| Action | File |
|--------|------|
| Create | `CRISPResso2/writers/__init__.py` (empty) |
| Create | `CRISPResso2/writers/vcf.py` (contents from CRISPRessoUtilities.py) |
| Create | `tests/unit_tests/test_writers/test_vcf/` directory |
| Modify | `CRISPResso2/CRISPRessoCORE.py` (update import + usage) |
| Move | `tests/.../test_CRISPRessoUtilities/*.py` to `tests/.../test_writers/test_vcf/` |
| Modify | Both test files (update imports, `utilities.` to `vcf.`) |
| Delete | `CRISPResso2/CRISPRessoUtilities.py` |
| Delete | `tests/unit_tests/test_CRISPRessoUtilities/` directory |
| Modify | `CLAUDE.md` (update documentation reference) |

## Implementation Steps

1. Create `CRISPResso2/writers/` directory with `__init__.py`
2. Copy `CRISPRessoUtilities.py` contents to `writers/vcf.py`
3. Update `CRISPRessoCORE.py` import and usage
4. Create test directory structure `tests/unit_tests/test_writers/test_vcf/`
5. Move test files and update their imports (`utilities.` to `vcf.`)
6. Delete old files (`CRISPRessoUtilities.py`, `test_CRISPRessoUtilities/`)
7. Update `CLAUDE.md` documentation
8. Run tests to verify nothing broke

## Design Decisions

**Why `writers/` instead of `utils/`:**
- `utils/` tends to become a catch-all for unrelated code
- `writers/` is specific to output file generation
- Other writers can be added here (BAM, JSON, etc.)

**Why not `CRISPRessoVCF.py` (flat structure):**
- There are other writers planned for this module
- `writers/` provides a clear organizational pattern for future additions

**No backwards compatibility:**
- Old import path (`from CRISPResso2 import CRISPRessoUtilities`) will break
- Acceptable because this is new code on a feature branch, not yet in production
