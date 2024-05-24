import xml.etree.ElementTree as ET
from pprint import pprint
def extract_data_from_xml(xml_file):
    # Parse the XML file
    tree = ET.parse(xml_file)
    root = tree.getroot()
    
    # List to store extracted data
    data = []
    
    # Iterate through each 'Hit' element
    for hit in root.findall('.//Hit'):
        # Extract 'Hit_def' and 'Hsp_evalue' elements
        hit_def = hit.find('Hit_def').text
        hsp_evalue = float(hit.find('.//Hsp_evalue').text)
        Hit_ascensionNumber=hit.find('./Hit_accession').text
        
        # Create a dictionary to store the data
        result_data = {'Hit_def': hit_def, 'Hsp_evalue': hsp_evalue, 'Ascension Number': Hit_ascensionNumber}
        
        # Append the dictionary to the list
        data.append(result_data)
    
    return data

# Example usage
xml_file_drought = 'resources/drought_resistance.xml'  # Replace 'example.xml' with the path to your XML file
extracted_data_drought = extract_data_from_xml(xml_file_drought)
print('DROUGHT RESISTANCE:')
pprint(extracted_data_drought[1:6])
xml_file_salinity = 'resources/salinity_tolerance.xml'  # Replace 'example.xml' with the path to your XML file
extracted_data_salinity = extract_data_from_xml(xml_file_salinity)
print('_'*200)
print('SALINITY TOLERANCE')
pprint(extracted_data_salinity[1:6])
print('_'*200)
xml_file_starch_quality = 'resources/starch_quality.xml'  # Replace 'example.xml' with the path to your XML file
extracted_data_starch_quality = extract_data_from_xml(xml_file_starch_quality)
print('STARCH QUALITY')
pprint(extracted_data_starch_quality[1:6])
print('_'*200)
