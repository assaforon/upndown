#  CDF used in many examples in an upcoming book on Up-and-Down by Flournoy and Oron

bookF = plogis(1:10,location=mloc,scale=slop)
# We assume 10 doses
round(bookF, 3)

