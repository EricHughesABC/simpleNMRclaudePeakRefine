# examplesPublic

Example exports that can be shared publicly. Only add datasets here that you have
permission to publish.

- `inputs/<name>.json`: the JSON exported from Mnova.
- `expected/<name>_refined.json`: the refined output it should give.

To create or update an expected output after checking it:

    simplenmr-claude-peak-refine examplesPublic/inputs/<name>.json examplesPublic/expected/<name>_refined.json

The `.mnova` documents themselves are never tracked.
