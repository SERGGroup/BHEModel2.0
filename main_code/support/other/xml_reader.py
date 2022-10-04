import xml.etree.ElementTree as ETree


def read_xml_file(file_path):

    with open(file_path, "rb") as xml_file:

        data = xml_file.read()

    root = ETree.fromstring(data)
    return root


def save_xml_file(root: ETree.Element, file_path):

    str_data = ETree.tostring(root)

    with open(file_path, "wb") as xml_file:

        xml_file.write(str_data)