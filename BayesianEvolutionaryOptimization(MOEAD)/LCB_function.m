function LCB = LCB_function(prefit,premse)
w = 2;
LCB = prefit-w*premse;
end