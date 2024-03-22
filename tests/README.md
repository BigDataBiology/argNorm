# Running the tests

## Test drug categorization
```
python -m pytest tests/test_arg_category.py
```

## Test normalizers
```
python -m pytest tests/test_normalizers.py
```

## Test CLI
```
bash tests/setup_cli_test.sh
python -m pytest tests/test_cli.py
```