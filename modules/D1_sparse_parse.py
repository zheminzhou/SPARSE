import sys, os, math, re


def report(level, fnames) :
    pathogens = {
        'Gardnerella ': '***',
        ' saprophyticus': '*** urinary tract infections',
        ' bovis': '** bovine pathogen',
        ' parasuis': '* swine pathogen',
        ' suis' : '** swine pathogen',
        'Leptospira ':'**** Leptospirosis',
        'Corynebacterium ulcerans':'***',
        'Mycobacterium ulcerans':'***',
        'Chlamydia ': '***',
        'enteric': '**',
        'Listeria monocytogenes':'* listeriosis',
        'Ureaplasma ':'*** urogenital diseases',
        'Ehrlichia ': '*** Ehrlichiosis',
        'Helicobacter ': '** Peptic ulcer',
        'Nocardia asteroides': '** Nocardiosis',
        'Francisella ': '*** Tularemia',
        'Chlamydophila ': '*** Psittacosis',
        'Mycobacterium avium' : '** animal pathogen',
        'tuberculosis':'*** TB',
        ' pestis':'**** plague',
        'Rickett': '*** Rocky mountain spotted fever',
        'Brucella ': '*** Brucellosis',
        ' pylori': '**** Peptic ulcer',
        'pertussis': '**** Whooping cough',
        'Bartonella ': '***',
        ' diphtheriae': '*** Diphtheria',
        'Shigella ': '**** Shigellosis',
        ' cholerae': '**** Cholera',
        'mallei': '** animal pathogen',
        ' leprae': '*** Leprosy ',
        ' africanum': '*** TB',
        ' gonorrhoeae': '** Gonorrhea',
        ' trachomatis': '*** Trachoma',
        ' canettii': '** TB',
        'pallidum': '**** syphilis',
        'Borrelia ': '**** Lyme disease/Relapsing fever',
        'Leishmania ': '**** Parasite',
        'Apicomplexa': '*** Parasite',
        ' pneumoniae': '** pneumonia',
        'Legionella pneumophila':'* legionellosis',
        'Salmonella ':'**** Salmonellosis/Typhoid',
        'Bacillus anthracis':'** Anthrax',
        'Campylobacter jejuni':'***',
        ' enterocolitica':'**',
        'Campylobacter coli':'***',
        'Staphylococcus aureus':'*',
        'Staphylococcus epidermidis':'*',
        'Pseudomonas aeruginosa':'*',
        'Moraxella ':'*',
        ' pyogenes':'*',
        'Escherichia ':'*',
        'Trichomonas vaginalis':'*** trichomoniasis',
        'Adenoviridae':'**** Virus',
        'Picornaviridae':'**** Virus',
        'Herpesviridae':'**** Virus',
        'Hepadnaviridae':'**** Virus',
        'Flaviviridae':'**** Virus',
        'Retroviridae':'**** Virus',
        'Orthomyxoviridae':'**** Virus',
        'Paramyxoviridae':'**** Virus',
        'Papovaviridae':'**** Virus',
        'Polyomavirus':'**** Virus',
        'Rhabdoviridae':'**** Virus',
        'Togaviridae':'**** Virus',
        'Filoviridae':'**** Virus',
        'Poxviridae':'**** Virus',
        'Neisseria meningitidis':'*** meningitis',
        'Burkholderia cenocepacia':'**',
        'Porphyromonas gingivalis':'* red_complex',
        'Tannerella forsythia':'* red_complex',
        'Treponema denticola':'* red_complex',
        'Aggregatibacter actinomycetemcomitans':'* periodontitis',
        'Capnocytophaga ':'* periodontitis',
        'Fusobacterium nucleatum':'* periodontitis',
        'Prevotella intermedia':'* periodontitis',
        'Entamoeba gingivalis':'* periodontitis',
        'Trichomonas tenax':'* periodontitis',
        'Fusobacterium polymorphum':'* periodontitis',
        'Fusobacterium necrophorum':'** Lemierre syndrome',
        'Streptobacillus moniliformis':'** Haverhill fever', 
        'Haemophilus influenzae':'* HACEK',
        'Haemophilus parainfluenzae':'* HACEK',
        'Haemophilus haemolyticus':'* HACEK',
        'Haemophilus parahaemolyticus':'* HACEK',
        'Aggregatibacter actinomycetemcomitans':'* HACEK',
        'Aggregatibacter segnis':'* HACEK',
        'Aggregatibacter aphrophilus':'* HACEK',
        'Aggregatibacter paraphrophilus':'* HACEK',
        'Cardiobacterium hominis':'* HACEK',
        'Cardiobacterium valvarum':'* HACEK',
        'Eikenella corrodens':'* HACEK',
        'Kingella kingae':'* HACEK',
        'Kingella denitrificans':'* HACEK',
        'Actinomyces israelii':'** Actinomycosis', 
        'Tropheryma whipplei':'* Whipple disease',
        'Arcanobacterium haemolyticum':'* Arcanobacterium haemolyticum infection',
        'Erysipelothrix rhusiopathiae':'* erysipeloid',
    }
    taxa = {'Unknown':['Dark Matter', 0.]}
    levels = level.split(',')
    data = {}
    fnames = [ fname if os.path.isfile(fname) else os.path.join(fname, 'profile.txt') for fname in fnames ]
    for fname in fnames :
        data[fname] = {}
        with open(fname) as fin :
            header = fin.readline().strip().split('\t')
            n_tot, n_found = int(header[1]), float(header[2])
            for line in fin :
                for level in levels :
                    if line.startswith(level) :
                        group, p1, p2, taxon = line.strip().split('\t')
                        taxon = taxon.rsplit(' (', 1)[0]
                        n_read = int(round(float(p2)/100. * n_found, 0))
                        if n_read > 0 :
                            n_tot -= n_read
                            data[fname][group] = n_read
                            if group not in taxa :
                                taxa[group] = [taxon, 0.]
                            taxa[group][1] += n_read
            data[fname]['Unknown'] = n_tot
    pathogens = sorted(pathogens.items(), key=lambda x:x[1])
    taxa_list = [t[0] for t in sorted(taxa.items(), key=lambda x:x[1][1], reverse=True)]
    print '#Group\t#Pathogenic\t{0}\t#Taxon\t#Pathogenic'.format('\t'.join(fnames))
    for group in taxa_list :
        label = 'non'
        taxon = taxa[group][0].rsplit('(', 1)[0]
        for p, t in pathogens :
            if re.findall(p, taxon) :
                label = t
                break
        print '{0}\t{3}\t{1}\t{2}'.format(group, '\t'.join([ str(data[fn].get(group, 0)) for fn in fnames ]), taxa[group][0], label)

if __name__ == '__main__' :
    report(sys.argv[1], sys.argv[2:])
