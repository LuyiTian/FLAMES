#parse config files
import json

def parse_json_config(json_file):
    with open(json_file) as f:
        decoded_dict = json.loads(f.read())
    # check if we have all the arguments:
    assert("transcript_counting" in decoded_dict)
    assert("isoform_parameters" in decoded_dict)
    assert("alignment_parameters" in decoded_dict)
    assert("global_parameters" in decoded_dict)
    assert("pipeline_parameters" in decoded_dict)
    assert(all(argument in decoded_dict["isoform_parameters"] for argument in ["MAX_DIST","MAX_TS_DIST","MAX_SPLICE_MATCH_DIST","Max_site_per_splice","Min_sup_cnt","Min_cnt_pct","Min_sup_pct","strand_specific","remove_incomp_reads"]))
    assert(all(argument in decoded_dict["alignment_parameters"] for argument in ["use_junctions"]))
    assert(all(argument in decoded_dict["global_parameters"] for argument in ["generate_raw_isoform"]))
    assert(all(argument in decoded_dict["pipeline_parameters"] for argument in ["do_genome_alignment","do_isoform_identification","do_read_realignment","do_transcript_quantification"]))
    assert(all(argument in decoded_dict["transcript_counting"] for argument in ["min_tr_coverage","min_read_coverage"]))
    # check all arguments within range
    assert(decoded_dict["isoform_parameters"]["MAX_DIST"]>0)
    assert(decoded_dict["isoform_parameters"]["MAX_TS_DIST"]>0)
    assert(decoded_dict["isoform_parameters"]["MAX_SPLICE_MATCH_DIST"]>0)
    assert(decoded_dict["isoform_parameters"]["Max_site_per_splice"]>0)
    assert(decoded_dict["isoform_parameters"]["Min_sup_cnt"]>0)
    assert(1>decoded_dict["isoform_parameters"]["Min_sup_pct"]>0)
    assert(decoded_dict["isoform_parameters"]["strand_specific"] in [-1,0,1])
    assert(decoded_dict["isoform_parameters"]["remove_incomp_reads"] >=0)

    assert(type(decoded_dict["alignment_parameters"]["use_junctions"]) == bool)

    assert(type(decoded_dict["global_parameters"]["generate_raw_isoform"]) == bool)
    assert(type(decoded_dict["global_parameters"]["has_UMI"]) == bool)
    return decoded_dict


def print_config(decoded_dict):
    print("Parameters in configuration file:")
    for cat1 in decoded_dict:
        if type(decoded_dict[cat1]) == dict:
            print(cat1)
            for cat2 in decoded_dict[cat1]:
                #print("\t",cat2,": ",decoded_dict[cat1][cat2],sep="")
                print "\t",cat2,":",decoded_dict[cat1][cat2]
        else:
            print cat1,":",decoded_dict[cat1]

if __name__ == "__main__":
    decoded_dict = parse_json_config("config_sclr_nanopore_default.json")
    print_config(decoded_dict)