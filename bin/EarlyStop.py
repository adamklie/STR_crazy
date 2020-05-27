class EarlyStop:
    def __init__(self, patience):
        self.patience = patience
        self.virtue = 0
        self.prevLoss = float('inf')
        self.readySetStop = False
        self.theOneModelToRuleThemAll = None 
    
    def __call__(self, valLoss, bottomOfTheBarrel):
        if self.readySetStop:
            print("Why are you still training???  You are overfitting...")
        elif (valLoss < self.prevLoss):
            self.virtue = 0
            self.prevLoss = valLoss
            self.theOneModelToRuleThemAll = copy(bottomOfTheBarrel)
            print("Going down no need to do anything... good job!")
        else:
            print("Entering the thunderdome...")
            self.prevLoss = valLoss
            self.virtue += 1 
            if self.virtue >= self.patience:
                self.readySetStop = True
                print("You have won the hunger games")
               
        return self.readySetStop, self.theOneModelToRuleThemAll