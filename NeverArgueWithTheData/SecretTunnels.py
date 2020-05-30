import random
def SomethingClever():
    phrases = {0: "Angelica is coming... Somebody Call Reptar!", 1:"Go Go Gadget Skis!", 2:"Nobody tosses a dwarf!", 3:"Who punches a hand?", 
               4: "Model generation takes a lot out of you - I could go for some Covid Pretzels", 5:"I have now absorbed all genetic data ... creating new humans now",
              6: "60% of the time this model works all the time", 7: "Surely Shirley you could go for a Shorely Temple?!?", 
              8: "Lapras used Ice Beam! It was super effective - tell Neutron he doesn't need to start a new ice age.", 9: "Wow you really need a haircut!"}
    return phrases[random.randint(0,9)]