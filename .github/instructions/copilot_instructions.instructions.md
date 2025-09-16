---
applyTo: '**/*.py' 
---
# Project coding standards

## General format
Python files should be formatted using black and isort. Documentation should be
compatible with Sphinx. Double quotes are preferred over single quotes.

## Docstrings
Function and Class docstrings should follow the following format:
```python
def function_name(param1: type, param2: type) -> type:
    """Brief description of the function.

    Longer description of the function, if necessary.

    Args:
        param1 (type): Description of param1.
        param2 (type): Description of param2.

    Raises:
        ExceptionType: Description of the exception raised, if any.

    Returns:
        type: Description of return value.
    """
```
## Comments
Since python black does not wrap comments, they should be wrapped to 88 characters as 
well.

## CLI arguments
CLI arguments should be parsed using argparse. All argparsing should be done in a 
separate function called `parse_args` that returns the dict of parsed arguments that
comes from `vars(parser.parse_args())`. 


---
applyTo: '**/*.ipynb'
---
Jupyter notebooks should be formatted using nbqa with black and isort. All other
formatting should follow the same rules as the Python files.