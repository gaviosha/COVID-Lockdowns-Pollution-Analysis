"""
Scrapes air quality data from European Environment Agecy site.
Run this script from ~data directory and adjust the global params below. POLLUTANT and POLLUTANT_CODE should be one of following tupples:

* {CO, 10}
* {NO, 38}
* {NO2, 8}
* {NOX as NO2, 9}
* {O3, 7}
* {PM10, 5}
* {SO2, 1}
"""

POLLUTANT = "O3"
POLLUTANT_CODE = 7
YEAR_FROM = 2018
YEAR_TO = 2021
HOW_SLOW = 0


import requests as r
import os
import re
import time
from string import Template
from tqdm import tqdm
from bs4 import BeautifulSoup


def make_data_dir_tree(pollutant, country_code, city_name, year_from, year_to):

    """
    Make all directories needed to store data.
    """

    home_dir = os.path.join(os.getcwd(), "raw")
    if not os.path.exists(home_dir):
        os.mkdir(home_dir)

    pollutant_dir = os.path.join(home_dir, pollutant)
    if not os.path.exists(pollutant_dir):
        os.mkdir(pollutant_dir)

    country_dir = os.path.join(pollutant_dir, country_code)
    if not os.path.exists(country_dir):
        os.mkdir(country_dir)

    city_dir = os.path.join(country_dir, city_name)
    if not os.path.exists(city_dir):
        os.mkdir(city_dir)

    for year in range(year_from, year_to + 1):
        year_dir = os.path.join(city_dir, str(year))
        if not os.path.exists(year_dir):
            os.mkdir(year_dir)

    return True


def get_data_urls(pollutant_code, country_code, city_name, year_from, year_to):

    """
    Get urls of all csv files containing pollution time series.
    """

    base_url = (
        "https://fme.discomap.eea.europa.eu/fmedatastreaming/AirQualityDownload/AQData_Extract.fmw?"
        "&CountryCode={CountryCode}"
        "&CityName={CityName}"
        "&Pollutant={Pollutant}"
        "&Year_from={Year_from}"
        "&Year_to={Year_to}"
        "&Station=&Samplingpoint=&Source=All&Output=HTML&UpdateDate=&TimeCoverage=Year"
    )

    url = base_url.format(
        CountryCode=country_code,
        CityName=city_name,
        Pollutant=pollutant_code,
        Year_from=year_from,
        Year_to=year_to,
    )

    page = r.get(url)
    soup = BeautifulSoup(page.text, "html.parser")

    urls_dict = dict()
    for year in range(year_from, year_to + 1):
        urls_dict[year] = []

    for ii in soup.find_all("a"):
        ii_url = ii.get("href")
        for year in range(year_from, year_to + 1):
            if "_" + str(year) + "_" in ii_url:
                urls_dict[year] += [ii_url]

    return urls_dict


def download_data_from_urls(
    urls_dict, pollutant, pollutant_code, country_code, city_name
):

    """
    Donload csv files one by one to correct folder.
    """

    home_dir = os.path.join(os.getcwd(), "raw")
    save_dir = os.path.join(home_dir, pollutant, country_code, city_name)

    str_start = str(pollutant_code) + "_"

    for year in urls_dict.keys():

        year_save_dir = os.path.join(save_dir, str(year))

        for href in urls_dict[year]:

            href_regex = re.escape(str_start) + r"(.*?)" + re.escape("_" + str(year))

            href_suffix = re.search(href_regex, href).group(1)

            href_file_name = (
                "_".join([pollutant, country_code, city_name, str(year), href_suffix])
                + ".csv"
            )

            href_file = r.get(href, allow_redirects=True)

            with open(os.path.join(year_save_dir, href_file_name), "wb") as f:
                f.write(href_file.content)

    return True


def download_data_wrapper(pollutant, pollutant_code, year_from, year_to, how_slow):

    home_dir = os.getcwd()
    all_countries_dir = os.path.join(
        home_dir, "countries_and_country_codes", "countries"
    )
    country_codes_dir = os.path.join(
        home_dir, "countries_and_country_codes", "country_codes"
    )

    all_countries = open(os.path.join(all_countries_dir, "countries.txt"))

    for country in tqdm(all_countries.readlines()):
        country_code = country[0:2]
        all_cities_loc = os.path.join(country_codes_dir, country_code + ".txt")

        if os.path.exists(all_cities_loc):
            all_cities = open(all_cities_loc)

            for city_name in all_cities.readlines():

                safe_city_name = city_name.strip("\n").replace("/", "")
                if safe_city_name == "":
                    pass

                else:
                    try:
                        time.sleep(how_slow)
                        make_data_dir_tree(
                            pollutant, country_code, safe_city_name, year_from, year_to
                        )
                        urls_dict = get_data_urls(
                            pollutant_code, country_code, city_name, year_from, year_to
                        )
                        download_data_from_urls(
                            urls_dict,
                            pollutant,
                            pollutant_code,
                            country_code,
                            safe_city_name,
                        )
                    except:
                        exception_str = ",".join(
                            [pollutant, country_code, safe_city_name]
                        )
                        with open("logger.txt", "a") as logger:
                            logger.write(exception_str + "\n")

            all_cities.close()

    all_countries.close()

    return True


if __name__ == "__main__":

    download_data_wrapper(POLLUTANT, POLLUTANT_CODE, YEAR_FROM, YEAR_TO, HOW_SLOW)
