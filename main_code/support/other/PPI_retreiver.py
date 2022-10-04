from main_code.constants import OTHER_FOLDER
from datetime import datetime
import requests
import warnings
import os


__DEFAULT_PPI = 239.79


def retrieve_PPI():

    """
        PPI is used to correct the well drilling cost according to the correlation provided by:

            Adams et al, “Estimating the Geothermal Electricity Generation Potential of Sedimentary Basins Using genGEO
            (the generalizable GEOthermal techno-economic simulator)", ChemRxiv Prepr., 2021.

            - if not explicitly set, PPI_current is AUTOMATICALLY RETRIEVED from:
                https://beta.bls.gov/dataViewer/view/timeseries/PCU2111--2111--

            - if the connection with the server is not possible, Jan 2022 data is used (PPI = 239.79)


        Once Updated PPI value is written on a specific file txt file in order to be retreived and checked by the
        researcher

    """

    __current_PPI = __get_stored_PPI()

    if __current_PPI is None:

        __current_PPI = __download_PPI()

        if not __current_PPI == __DEFAULT_PPI:

            __store_updated_PPI(__current_PPI)

    return __current_PPI

def __get_stored_PPI():

    currentMonth = datetime.now().month
    currentYear = datetime.now().year
    overall_month = currentYear * 12 + currentMonth

    try:

        with open(os.path.join(OTHER_FOLDER, "PPI_value.txt"), "r") as file:

            line = file.readline().strip("\n")

        split_line = line.split(" - ")
        split_date = split_line[0].split("/")
        PPI = float(split_line[1])

        delta_month = overall_month - (int(split_date[1]) * 12 + int(split_date[0]))

        if delta_month < 3:

            return PPI

        else:

            return None

    except:

        return None

def __download_PPI():

    __current_PPI = __DEFAULT_PPI

    r = requests.get('https://api.bls.gov/publicAPI/v2/timeseries/data/PCU2111--2111--', params={'latest': 'true'})

    if r.status_code == 200:

        try:

            content = r.json()
            __current_PPI = float(content["Results"]["series"][0]['data'][0]['value'])

        except:

            warnings.warn(

                "Impossible to retrieve LATEST PPI value from BLS database,\nDefault PPI used "
                "instead:\ncurrent_PPI = {}".format(__DEFAULT_PPI)

            )

    return __current_PPI

def __store_updated_PPI(PPI):

    currentMonth = datetime.now().month
    currentYear = datetime.now().year

    with open(os.path.join(OTHER_FOLDER, "PPI_value.txt"), "w") as file:

        file.write("{}/{} - {}".format(int(currentMonth), int(currentYear), PPI))
