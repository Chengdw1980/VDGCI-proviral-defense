import re
DIRECT={'Cas','McrBC','GmrSD_RM_Type_IV','Other_Type_IV','EcoKMcrA','ScoMcrA','MspJI','DpnI','PvuRts1I','BREX','DISARM','Dnd','SspBCDE','pAgo','Wadjet','DdmDE','SNIPE','TagI'}
ABITA={'Lit','PrrC','RloC','MazEF','RosmerTA','PfiAT','DarTG','RnlAB','SanaTA','PsyrTA','TgvAB','PhiR47_GevAT','Mok_Hok_Sok','CapRel','Retron','CBASS','Pycsar','Thoeris','Avs','CARD_NLR','NLR','Rst_TIR-NLR','Paris','ToxIN','AbiEii'}
def mech(t):
    if t.startswith('RM_') or t in DIRECT: return 'Direct nucleic-acid targeting'
    if t.startswith('Abi') or t in ABITA: return 'Abortive infection / toxin-antitoxin'
    return 'Other or mechanism unresolved'
