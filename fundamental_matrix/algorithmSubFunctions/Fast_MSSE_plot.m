function Fast_MSSE_plot( fiterrs, Lambda, m_p)

    numpts = length(fiterrs);
    
    R2= sort(fiterrs.^2);
    msse_crit = ((1 : numpts)'-m_p).*R2./cumsum(R2);
    msse_crit(1:m_p) = 1;
        
    plot(msse_crit.^0.5,'-*'), hold on
    plot([1 numpts],[Lambda Lambda],'g')
    hold off
end