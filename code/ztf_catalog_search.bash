#catalog search script with table upload
#can edit this for any catalog (IE COSMOS)

curl -F filename=@merian_pos_ipac.dat \
-F catalog=ztf_objects_dr22 \
-F spatial=Upload \
-F uradius=2 \
-F outfmt=1 \
"https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query" \
-o merian_ztf_xmatch_output.tbl 