import random

def peek(S):
	x = S.pop()
	S.add(x)
	return x

def weighted_choice(choices):
   total = sum(w for c, w in choices)
   r = random.uniform(0, total)
   upto = 0
   for c, w in choices:
      if upto + w >= r:
         return c
      upto += w
   print choices
   assert False, "Shouldn't get here"

def reverse_string(s):
   return s[::-1]
