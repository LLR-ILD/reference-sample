from pysteer import Pysteer

# -----------------------------------------------------------------------------
# Some basic test/use examples for pysteer.
if __name__ == "__main__":
    cpd = {"InitializeDD4hep": {"DD4hepXMLFile": {"value": "dump.shit"},
                                        "EncodingStringParameter": "c"
                                        },
        "non-ex-pr": {"this": "should never be read"}
    }
    steerer = Pysteer(change_parameter_defaults=cpd,
        set_parameter_value={"EncodingStringParameterName": "stringper"})
    steerer.add("TauConesProcessor", {"IsolationConeAngle": dict(value=".3")})
    steerer.add("TauConesProcessor")
    steerer.add("IsolatedLeptonTaggingProcessor")
    steerer.write(xml_name="steer.xml")
    print(steerer)
    print(steerer.marlin_global)
