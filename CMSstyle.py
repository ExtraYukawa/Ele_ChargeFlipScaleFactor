from ROOT import *

def SetStyle(gPad, era):

        if(era == '2016postapv'):
          era = '2016' # To have a more beautiful plot
	latex = TLatex();
	latex.SetNDC()
	l = gPad.GetLeftMargin();
	t = gPad.GetTopMargin();
	r = gPad.GetRightMargin();
	b = gPad.GetBottomMargin();
	#CMS text
	cmsText = "CMS";
	cmsTextFont = 60;
	cmsTextSize = 0.6;
	cmsTextOffset = 0.1;
	
	relPosX = 0.12;
	relPosY = 0.035;
	
	# extra 
	extraText = "  Preliminary "+era;
	#extraText = "";
	extraOverCmsTextSize = 0.76;
	extraTextFont = 52
	
	lumiText = ""
        if(era=='2017'):
          lumiText = "41.5 fb^{-1} (13 TeV)"
        elif(era == '2018'):
          lumiText = "59.8 fb^{-1} (13 TeV)"
        elif(era == '2016apv'):
          lumiText = "19.5 fb^{-1} (13 TeV)"
        elif(era == '2016'):
          lumiText = "16.8 fb^{-1} (13 TeV)"
	lumiTextSize = 0.5;
	lumiTextOffset = 0.2;
	relExtraDY = 1.2;
	
	
	latex.SetTextAngle(0);
	latex.SetTextColor(kBlack);
	extraTextSize = extraOverCmsTextSize * cmsTextSize;
	latex.SetTextFont(42);

	latex.SetTextFont(42);
	latex.SetTextAlign(31);
	latex.SetTextSize(lumiTextSize * t);
	latex.DrawLatex(1 - r, 1 - t + lumiTextOffset * t, lumiText);
	
	latex.SetTextFont(cmsTextFont);
	latex.SetTextAlign(11);
	latex.SetTextSize(cmsTextSize * t);
	latex.DrawLatex(l, 1 - t + lumiTextOffset * t, cmsText);
	
	posX_ = 0
	posX_ = l + relPosX * (1 - l - r);
	posY_ = 1 - t - relPosY * (1 - t - b);
	posX_ = l + relPosX * (1 - l - r);
	posY_ = 1 - t + lumiTextOffset * t;
	alignX_ = 1;
	alignY_ = 1;
	align_ = 10 * alignX_ + alignY_;
	latex.SetTextFont(extraTextFont);
	latex.SetTextSize(extraTextSize * t);
	latex.SetTextAlign(align_);
	latex.DrawLatex(posX_, posY_, extraText);
	latex.SetTextAlign(31);
	latex.SetTextSize(lumiTextSize * t);
	return gPad
