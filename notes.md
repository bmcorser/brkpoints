# Good
- The code might work?

# Bad
- No packaging integration
- Can't run tests without 10GB of data?
- Running tests as described in docs fails with `ImportError`
- Code doesn't compile `AttributeError: _Environ instance has no __call__ method` cf7140f
- Code doesn't pass PEP8
- [Illegible syntax](https://github.com/DeskGen/brkpoints/pull/5#discussion_r24162577)
- Single character variable names in nested loops aren't readable even if
  they're not in a convoluted list comp
- "Predictably" formatted strings used instead of
  [data](https://docs.python.org/2/library/stdtypes.html#dict)
  [structures](https://docs.python.org/2/library/stdtypes.html#typesseq)!?
- Reticent to assign? 438f80a
