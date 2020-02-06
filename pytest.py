import pysteer
from pysteer import marlin_processors_dict, update_registered

global_default_dict = {
    '"LCIOInputFiles"': dict(value=(
        "\n            ".join(["",
            #"/home/kunath/ILD/Data_SM/higgs_ffh/rv01-16-p10_250.sv01-14-01-p00.mILD_o1_v05.E250-TDR_ws.I106485.Pqqh.eL.pR-00001-DST.slcio",
            "/home/kunath/ILD/Data_SM/higgs_ffh/rv01-16-p10_250.sv01-14-01-p00.mILD_o1_v05.E250-TDR_ws.I106482.Pe3e3h.eR.pL-00001-DST.slcio",
    ])+"\n   ")),
    '"MaxRecordNumber"': dict(value=1000),
    '"SkipNEvents"': dict(value=844),
    '"Verbosity"': dict(value="MESSAGE0",
        options = "DEBUG0-4,MESSAGE0-4,WARNING0-4,ERROR0-4,SILENT"),
}

# -----------------------------------------------------------------------------
def xml_parameters(default_dict, changes_dict={}):
    param_strings = []
    for parameter_name, parameter_dict in sorted(default_dict.items()):
        if (parameter_name in changes_dict
        and changes_dict[parameter_name]["value"] != parameter_dict["value"]):
            parameter_dict = changes_dict[parameter_name]
            indent = "    "
        else: # Indent the default parameters differently.
            indent = "        "
        param_line = indent + "<parameter name={}".format(parameter_name)
        for attr in parameter_dict:
            if attr == "value": continue
            param_line += " {tag}={value}".format(
                tag=attr, value=parameter_dict[attr])
        param_line += "> {} </parameter>".format(parameter_dict["value"])
        param_strings.append(param_line)
    return "\n".join(param_strings)





def marlin_xml(project_defaults, global_dict={}, processors=[], xml_name="steering.xml"):
    """Produce an xml file that can be fed into Marlin.

    The global dict can be filled with the global values that we want to alter
    wrt the defaults.
    Processors expects the ordered list of processor tuples. The first entry in
    the tuple should be the processor name, the second (and last) entry should
    be a dict with those parameters that we want to alter wrt the defaults (can
    be an empty dict).
    """
    header_comment = ("<!-- This steering file was produced via a python "
      "script.-->\n\n")
    xml_pieces = []
    xml_pieces.append(header_comment)
    xml_pieces.append("<marlin>\n\n  <execute>\n")
    called_processor_types = list(list(zip(*processors))[0])
    called_processor_names = []
    for i, processor_type in enumerate(called_processor_types):
        name = processor_type[:-1]+'{:03}"'.format(i)
        called_processor_names.append(name)
        xml_pieces.append("    <processor name={}/>\n".format(name))
    xml_pieces.append("  </execute>\n\n  <global>\n")
    xml_pieces.append(xml_parameters(
        default_dict=global_default_dict,
        changes_dict={}
    ))
    xml_pieces.append("\n  </global>\n\n")

    for name, proc in zip(called_processor_names, processors):
        proc_type, proc_dict = proc
        header = "  <processor name={} type={}>\n".format(name, proc_type)
        xml_pieces.append(header)
        xml_pieces.append(xml_parameters(
            default_dict=project_defaults[proc_type],
            changes_dict=proc_dict
        ))
        xml_pieces.append("\n  </processor>\n\n")
    xml_pieces.append("</marlin>")
    xml_content = "".join(xml_pieces)
    with open(xml_name, "w") as write_file:
        write_file.write(xml_content)







if __name__ == "__main__":
    steerer = pysteer.Pysteer()
    steerer.update_processors()
    print(steerer.ilcsoft_processors)
    #
    changes_dict={
        '"IsolationConeAngle"': dict(value=".3"),
        '"IsolationEnergy"': dict(value="2.0"),
        '"SearchConeAngle"': dict(value="0.125"),
        }
    marlin_xml(
        project_defaults=
        processors=[
        ('"TauConesProcessor"', changes_dict),
        #('"TauConesProcessor"', changes_dict)
    ])
