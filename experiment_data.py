import numpy as np
import json
import os

from glob import glob
from h5py import File


class ExperimentData:

    def __init__(self, hdf_path):

        hdf = File(hdf_path, "r")

        self.policies = list(hdf.keys())
        self.award_amounts = list(
            hdf[self.policies[0]].keys()
        )
        self.pubneg_rates = list(
            hdf[self.policies[0]][self.award_amounts[0]].keys()
        )
        self.fpdrs = list(
            hdf[self.policies[0]][
                    self.award_amounts[2]][
                        self.pubneg_rates[0]].keys()
        )

        self.hdf = hdf

    @classmethod
    def from_jsons(cls, jsons_dir, hdf_path):
        """
        Create new instance from a directory of jsons (raw output from
        scimod-agency). Write the new HDF to hdf_path. HDF stays open unless
        closed explicitly. Using Python destructor to close in case user
        forgets.
        """
        hdf = File(hdf_path, "w")

        json_files = glob(os.path.join(jsons_dir, "*"))
        n_json = len(json_files)
        for idx, jf in enumerate(json_files):
            print("processing JSON {}/{}".format(idx + 1, n_json))
            j = json.load(open(jf, "r"))
            md = j["metadata"]
            policy = md["policy"]

            parameters = md["parameters"]
            award_amount = parameters["awardAmount"]
            pubneg_rate = parameters["publishNegativeResultRate"]
            fpdr = parameters["falsePositiveDetectionRate"]

            group = "/{}/{}/{:.2f}/{:.2f}".format(
                policy, award_amount, pubneg_rate, fpdr
            )
            try:
                hdf.create_group(group)
                for measure in [
                            'falseDiscoveryRate',
                            'falsePositiveRate',
                            'meanPublications',
                            'sumPublications',
                            'medianPublications',
                            'sumFunds',
                            'medianFunds',
                            'meanFunds'
                        ]:
                    try:
                        hdf[group].create_dataset(
                            measure, data=np.array(j[measure], dtype='float'),
                            compression="gzip", compression_opts=9
                        )
                    except KeyError:
                        pass

                if len(j['agentFPRs']) > 0:
                    hdf[group].create_dataset(
                        'agentFPRs', data=np.array(j['agentFPRs'], dtype='float'),
                        compression="gzip", compression_opts=9
                    )


            except ValueError:
                pass

        hdf.close()

        return cls(hdf_path)

    def __getitem__(self, key):

        if type(key) is tuple:
            path = ("/{}"*len(key)).format(*key)
        elif type(key) is str:
            path = "/" + key

        return self.hdf[path]

    def __del__(self):
        self.hdf.close()
