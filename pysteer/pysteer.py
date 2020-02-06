################################################################################
#
#  Class definition file for Pysteer.
#
################################################################################
#from .extract_from_cpp import registered_parameter_dict
from .write_processor_parameters import update_registered

class Pysteer(object):
    """Interface for the creation of Marlin steering files.

    For now, just builds a dictionary file with info on the registered
    parameters in the Marlin processors that are linked with the class.

    :param confirm_ilcsoft_defaults [bool]: If false, the dict of external
        processor parameter information is only updated if a processor is
        requested that is not yet in the dict.
    : param ilcsoft_path [str]: Location from which Marlin will be run.
    : param ilcsoft_processors [list(str)|None]: Defines which (external) Marlin
        processors this class should know about. If None is chosen, all
        processors available in the chosen Marlin version will be available, for
        the price of reduced clarity.
    """
    def __init__(
        self,
        confirm_ilcsoft_defaults=False,
        ilcsoft_path="/cvmfs/ilc.desy.de/sw/x86_64_gcc49_sl6/v02-00-02",
        ilcsoft_processors=[
            "InitializeDD4hep",
            "IsolatedLeptonTaggingProcessor",
        ],
    ):
        self.confirm_ilcsoft_defaults=confirm_ilcsoft_defaults
        self.ilcsoft_path=ilcsoft_path
        self.ilcsoft_processors=ilcsoft_processors

        self.update_processors()

    def update_processors(self):
        """Get the .json files in `pysteer` that store information on the
        default values of project and external Marlin processors up-to-date.

        If self.confirm_ilcsoft_defaults=False, the dict of external
        processor parameter information is only updated if a processor is
        requested that is not yet in the dict.
        """
        update_registered(
            confirm_ilcsoft_defaults=self.confirm_ilcsoft_defaults,
            ilcsoft_path=self.ilcsoft_path,
            load_only=self.ilcsoft_processors,
        )

    def full_update(self):
        """Rebuild the processor default parameter dictionaries.

        Equivalent to self.update_processors() if
        self.confirm_ilcsoft_defaults=True.
        """
        update_registered(
            confirm_ilcsoft_defaults=True,
            ilcsoft_path=self.ilcsoft_path,
            load_only=self.ilcsoft_processors,
        )



