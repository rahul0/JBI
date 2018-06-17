
# Abstract
In clinical practice, the success of clinical trial designs, association analysis, documentation
of mandatory discharge summaries, and other tasks depends on answering
questions such as “Has this person ever been treated for breast cancer?” This thesis
considers the limited success of traditional machine learning approaches with respect
to addressing these problems, and proposes new methods for answering such questions.
Our solution is to develop methods that automatically highlight key textual passages
to provide help in answering the questions. The advantages of this approach lie
in the fact that it does not require experts to go through the whole dataset of electronic
medical records. Without such an annotation, an answer to “Has this person
ever been treated for breast cancer” would otherwise only be obtained when reading
manually through the whole electronic text of the medical record.



## Shorter version: http://cmj4.web.rice.edu/TACL.pdf
## Longer version aka master thesis: https://scholarship.rice.edu/handle/1911/77189

Nomenclature are very closed to as discussed in paper/thesis.

# JBI
JBI paper and few processing script.

### models/
 ├── 1a  :WLR-no-expert MD_A

 ├── 1b  :WLR-no-expert MD_B

 ├── 2a  :WLR-expert MD_A

 ├── 2b  :WLR-expert MD_B

 └── 3a  :NOT SUPPORTED, This model has been discussed in thesis but not in paper.


### scripts/

 ├── bashrc_sample  : compilation varible to compile with minuit2 and GSL has been defined here.

 ├── compile.sh : To compile model, have compiler optimization flags etc.

 ├── compile_underlining.sh: To compile underlining model.

 └── makeHTML.sh : I used to convert annotated text in html format to get better visual interpretation.

 Please feel free to reach if you would like certain explaination or would like to incorporate code in software.
