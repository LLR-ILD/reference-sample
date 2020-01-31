import re
def remove_c_comments(text):
    """ Remove c-style comments.

    text: Blob of text with comments (can include newlines).
    returns: Text with comments removed.
    Found at:
        https://gist.github.com/ChunMinChang/88bfa5842396c1fbbc5b.
    """
    replacer = lambda match: (
        " " if match.group(0).startswith('/') else
        match.group(0)
    )
    pattern = re.compile(
        r'//.*?$|/\*.*?\*/|\'(?:\\.|[^\\\'])*\'|"(?:\\.|[^\\"])*"',
        re.DOTALL | re.MULTILINE
    )
    return re.sub(pattern, replacer, text)

def expression_end_position(code, char=";"):
    """Returns position of the first semicolon outside of a string.

    Assumes that there are no comments in the string!
    """
    semicolon_position = code.find(char)
    without_strings = code[:semicolon_position]
    semicolon_is_accepted = False
    while not semicolon_is_accepted:
        while ((without_strings.find("'") != -1)
        or (without_strings.find("\"") != -1)):
            if without_strings.find("'") * without_strings.find("\"") <= 0:
                # Only one of the finds was successful.
                string_starts_at = max(
                    without_strings.find("'"), without_strings.find("\""))
            else:
                string_starts_at = min(
                    without_strings.find("'"), without_strings.find("\""))
            string_length = without_strings[string_starts_at+1:].find(
                without_strings[string_starts_at]) + 2
            if string_length == 1:
                # This means that the char was inside a comment.
                semicolon_position += code[semicolon_position+1:].find(char) +1
                break
            without_strings = (without_strings[:string_starts_at]
                + without_strings[string_starts_at+string_length:])
        semicolon_is_accepted = True
    return semicolon_position

def c_function_parameter_list(code, function_name):
    """Returns the list of function parameters of the first occurance of the
    funtion as well as the end position of this first occurance.

    Returns None object and position -1 if the function never occures.
    Assumes that there are no comments in the string!
    """
    function_begin = code.find(function_name)
    if function_begin == -1:
        return None, -1
    while code[function_begin+len(function_name)] not in " (":
        try_new_find = code[function_begin+1:].find(function_name)
        if try_new_find == -1:
            return None, -1
        function_begin += try_new_find + 1

    from_func_till_end = code[function_begin:]
    func_end_position = expression_end_position(from_func_till_end)
    function = from_func_till_end[:func_end_position]
    open_br = function.find("(")
    close_br = function.rfind(")") # rfind: Last occurance.
    if open_br == -1 or close_br == -1:
        return None, -1
    within_brackets = function[open_br+1:close_br]
    ## Below line did not ignore commas inside of comments.
    #parameter_list = within_brackets.split(",")
    parameter_list = []
    comma_position = -1
    length_of_parameter = -1
    while length_of_parameter != 0:
        length_of_parameter = expression_end_position(
            within_brackets[comma_position+1:], char=",") + 1
        parameter_list.append(within_brackets[
            comma_position+1:comma_position+length_of_parameter])
        comma_position += length_of_parameter
    # The last list item is empty. Fill it with the content after the final
    parameter_list[-1] = within_brackets[comma_position+1:]
    for i, param in enumerate(parameter_list):
        parameter_list[i] = param.rstrip().lstrip()
    # Remove the last item if it is empty (trailing comma in code).
    if parameter_list[-1] == "":
        parameter_list.pop()
    return parameter_list, function_begin+func_end_position

def c_function_parameter_all_occurances(code, function_name):
    """Returns a  list of all the parameter lists for a functions occurance in
    the code. (Empty list for no occurance).

    Assumes that there are no comments in the string!
    """
    parameter_list_list = []
    remaining_code = code
    while True:
        parameter_list, start = c_function_parameter_list(
            remaining_code, function_name)
        remaining_code = remaining_code[start:]
        if parameter_list:
            parameter_list_list.append((parameter_list))
        else:
            break
    return parameter_list_list

def cpp_value_parser(value):
    """Strip some common cpp decorators. E.g.
    5.0f -> 5.0, std::string("abc") -> "abc"
    """
    if value[-1] == "f":
        return value[:-1]
    if 'std::string(' in value:
        value =  value.replace("std::string(", "")[1:-2]
    return value

def marlin_register_dict(code):
    """Create a dict of the registered parameters and collections of a Marlin
    processor (ILCSoft).

    Assumes that there are no comments in the string!
    """
    # For other things we might want to specify. E.g. lcioInType, type.
    description_dict = {}
    param_dict = {}
    processor_parameters = c_function_parameter_all_occurances(
        code, function_name="registerProcessorParameter")
    for parameter in processor_parameters:
        assert len(parameter) == 4 # [name, description, cpp_var_name, value]
        name = parameter[0]
        description = parameter[1]
        value = cpp_value_parser(parameter[3])
        description_dict[name] = description
        param_dict[name] = dict(value=value)

    processor_parameters = c_function_parameter_all_occurances(
        code, function_name="registerInputCollection")
    for parameter in processor_parameters:
        # [lcioInType, name, description, cpp_var_name, value]
        assert len(parameter) == 5
        lcioInType = parameter[0]
        name = parameter[1]
        description = parameter[2]
        value = cpp_value_parser(parameter[4])
        description_dict[name] = description
        param_dict[name] = dict(value=value, lcioInType=lcioInType)

    processor_parameters = c_function_parameter_all_occurances(
        code, function_name="registerOutputCollection")
    for parameter in processor_parameters:
        # [lcioOutType, name, description, cpp_var_name, value]
        assert len(parameter) == 5
        lcioOutType = parameter[0]
        name = parameter[1]
        description = parameter[2]
        value = cpp_value_parser(parameter[4])
        description_dict[name] = description
        param_dict[name] = dict(value=value, lcioOutType=lcioOutType)

    return  param_dict, description_dict

def marlin_processors_dict(look_at_list):
    """
    """
    param_dict = {}
    descriptions_dict = {}
    for cc_file_name in look_at_list:
        with open(cc_file_name, "r") as file:
            comment_free_cc = remove_c_comments(file.read())

        processor_name = c_function_parameter_list(comment_free_cc,
                                                   "marlin::Processor")[0][0]
        param_1proc, descriptions_1proc = marlin_register_dict(comment_free_cc)
        param_dict[processor_name] = param_1proc
        descriptions_dict[processor_name] = descriptions_1proc
    return param_dict

# -----------------------------------------------------------------------------
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
cc = "/home/kunath/reference-sample/processors/tau_cones/src/tau_cones_processor.cc"
project_processors = marlin_processors_dict([cc])

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





def marlin_xml(global_dict={}, processors=[], xml_name="steering.xml"):
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
            default_dict=project_processors[proc_type],
            changes_dict=proc_dict
        ))
        xml_pieces.append("\n  </processor>\n\n")
    xml_pieces.append("</marlin>")
    xml_content = "".join(xml_pieces)
    with open(xml_name, "w") as write_file:
        write_file.write(xml_content)







if __name__ == "__main__":
    parameter = "<parameter >"
    #print(parameter.format("Reader", "code"))

    changes_dict={
        '"IsolationConeAngle"': dict(value=".3"),
        '"IsolationEnergy"': dict(value="2.0"),
        '"SearchConeAngle"': dict(value="0.125"),
        }

    #for processor_name in project_processors:
    #    print(xml_parameters(
    #        default_dict=project_processors[processor_name],
    #        changes_dict=changes_dict
    #    ))
    marlin_xml(processors=[
        ('"TauConesProcessor"', changes_dict),
        #('"TauConesProcessor"', changes_dict)
    ])

