import enum

class Shorthands(enum.Enum):
    HALOGEN = '[F,Cl,Br,I]'
    AROMATIC = 'a'
    ALIPHATIC = '[C,N,O,S,P]'
    ANY = '*'  # does not match hydrogen
    DELETE = ''
    ALIPHATIC_DONOR = '[N,O,S;!H0;v3,v4&+1]'
    AROMATIC_DONOR = '[n,o,s;!H0;+0]'
    ALIPHATIC_ACCEPTOR = '[N,O,S;H0;v2;!$(*-*=[O,N,P,S])]'
    AROMATIC_ACCEPTOR = '[n,o,s;H0;+0;!$(*-*=[O,N,P,S])]'
