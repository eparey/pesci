"""
Global definition for tqdm progress bar format.
"""


def supports_color():
    """
    Returns True if the running system's terminal supports color, and False
    otherwise.
    """
    plat = sys.platform
    supported_platform = plat != 'Pocket PC' and (plat != 'win32' or
                                                  'ANSICON' in os.environ)
    # isatty is not always implemented, #6223.
    is_a_tty = hasattr(sys.stdout, 'isatty') and sys.stdout.isatty()
    return supported_platform and is_a_tty

if supports_color:
    BAR_FORMAT = '{percentage:3.0f}% |{bar:50}| task \x1B[1;32m{unit}'
else:
    BAR_FORMAT = '{percentage:3.0f}% |{bar:50}| task {unit}'
