################################################################################
#
#  A mini-class for conveniently changing the <global/> section of a Marlin
#  steering file build with the pysteer class.
#
################################################################################
class MarlinGlobal(object):
    """Object that stores the information for the <global/> section.
    Most likely LCIOInputFiles will have to be changed.

    : param MaxRecordNumber (int|str): If 0, all available events are used.
        Else, stop after MaxRecordNumber events. Default=0
    : param SkipNEvents (int|str):  Default=0
    : param Verbosity (str): Valid options: DEBUG0-4, MESSAGE0-4, WARNING0-4,
        ERROR0-4,SILENT. Default=MESSAGE.
    : param LCIOInputFiles (lis[str]): Location of the .slcio files to read.
        It is very unlikely that the default valu is what you want.
    """
    def __init__(
        self,
        MaxRecordNumber=0,
        SkipNEvents=0,
        Verbosity="MESSAGE",
        LCIOInputFiles="\n          ".join(["",
                "/home/kunath/ILD/Data_SM/higgs_ffh/rv01-16-p10_250.sv01-14-01-p00.mILD_o1_v05.E250-TDR_ws.I106485.Pqqh.eL.pR-00001-DST.slcio",
                "/home/kunath/ILD/Data_SM/higgs_ffh/rv01-16-p10_250.sv01-14-01-p00.mILD_o1_v05.E250-TDR_ws.I106482.Pe3e3h.eR.pL-00001-DST.slcio",
        ])+"\n     ",
    ):
        self.LCIOInputFiles = LCIOInputFiles
        self.MaxRecordNumber = MaxRecordNumber
        self.SkipNEvents = SkipNEvents
        self.Verbosity = Verbosity

    def as_dict(self):
        """Build a python dict version of the global information.
        """
        global_dict = dict(
            LCIOInputFiles=dict(value=self.LCIOInputFiles),
            MaxRecordNumber=dict(value=self.MaxRecordNumber),
            SkipNEvents=dict(value=self.SkipNEvents),
            Verbosity=dict(value=self.Verbosity,
                options="DEBUG0-4,MESSAGE0-4,WARNING0-4,ERROR0-4,SILENT"),
        )
        return global_dict

    def __str__(self):
        gl_print = ["A MarlinGlobal object with values:"]
        for param, value in self.__dict__.items():
            gl_print.append("{}: {}".format(param, value))
        return "\n  ".join(gl_print)



