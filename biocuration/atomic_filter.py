from biocuration import uniprot as up
from biocuration.uniprot.atomic import APile


p = APile()
with open('/home/klemens/Downloads/allnew.txl', 'r', encoding='ascii') as infile:
    for entry in up.parse_txt_compatible(infile):
        p.consume(entry)
    print(p.size())