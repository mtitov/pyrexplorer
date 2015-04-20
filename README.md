pyReXplorer
=====

pyReXplorer (Relations eXplorer) is a python package aimed to explore relations 
in analyzed data (patterns discovery, i.e. data mining).

Sequential Pattern Mining

- Algorithm SPADE (Sequential PAttern Discovery using Equivalence classes). This 
implementation is based on paper “SPADE: An Efficient Algorithm for Mining 
Frequent Sequences” by Mohammed J. Zaki and was encouraged by Shahin Saneinejad 
(github.com/shahin) with “sequenceminer”.



Examples
=====

SPADE

    from pyrexplorer.spade import SPADE

    input_data = {
        'sequences': {
            1: {2: ('A', 'B', 'C'), 4: ('B', 'C')},
            2: {3: ('B', 'C')}
        },
        'minimum_support': 2
    }

    spade = SPADE()
    spade.set(**input_data)
    for element in spade.execute(sort=False, k=None, top_number=2):
        print "k={0:<8}supp={1:<10}seq={2}".format(len(element),
                                                   element.support,
                                                   element.sequence)
