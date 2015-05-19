
files='pE_169 pE_216 pE_356 pE_368 pE_499 pJ_169 pJ_216 pJ_356 pJ_368 pJ_499 p_128 p_141 p_149 p_191 p_20 p_21 p_210 p_214 p_22 p_23 p_25 p_26 p_266 p_29 p_30 p_303 p_31 p_34 p_35 p_36 p_37 p_38 p_39 p_40 p_41 p_42 p_422 p_43 p_44 p_45 p_46 p_47 p_477 p_48 p_49 p_50 p_51 p_53 p_55 p_56 p_96 p_97 p_98'.split()
#files='pE_169 pE_216'.split()
openingTag='<:t:>'
closingTag='<:/t:>'

def extract_phrase_between_tags( Record):
         """ find value between tags, if that tag is not present returns None"""

         if( not (Record.find(openingTag)!=-1 and Record.find(closingTag)!=-1) ):
             return None
         temp1= Record.find(openingTag)
         temp2= Record.find(closingTag)

         return Record[temp1+len(openingTag):temp2]


def extract_tagged_token_count(sentence, filename ):
   """extract count of tagged tokens """
   # discard punctuations
   punctuations=set([',','.', openingTag, closingTag])
   if sentence=='':
      return


   File2String=sentence.replace( openingTag.upper(),openingTag ).replace( closingTag.upper(), closingTag)
   startIndex=File2String.find(openingTag)

   #This means document do not have the tag. I dont think that is poosible, so test this hyposthesis.
   assert startIndex!=-1
   token_count=0

   while startIndex!=-1:

        #find where tag ends
        end_index=File2String.find(closingTag)
        #assert end_index!=-1, File2String
        #grab every thing between tags"""
        fne_phrase=extract_phrase_between_tags(' '.join(File2String.split()).lower().strip())

        # incement the length
        if  fne_phrase:
            # clean if any  opening or closing tag present
            fne_phrase=fne_phrase.replace(openingTag,'').replace(closingTag,'')
            tokens=fne_phrase.split(' ')
            token_count+= len([ token for token in tokens if token.strip() not in punctuations and len(token)  ])
            #our new string is truncated one"""
            File2String=File2String[end_index+ len(closingTag):]
            #go to this new position where fne starts"""
            startIndex=File2String.find(openingTag)
        else:
            print 'Error in file:'.format(filename)


   assert token_count>0
   return token_count



if __name__=='__main__':

    for  f in files:
        content=' '.join(open(f).readlines())
        print'{}:{} '.format(f,extract_tagged_token_count(content, f))

else:
    pass


