function [dfdu_full_abs] = dfdu_full_abs(x,p)

  dfdu_full_abs(1,1)=0;
  dfdu_full_abs(2,1)=0;
  dfdu_full_abs(3,1)=(p.a_len^2*p.m_leg + p.l^2*p.m_leg + p.l^2*p.mh - p.b_len*p.l*p.m_leg*...
         cos(x(1) - x(2)))/(p.b_len^2*p.m_leg*(p.a_len^2*p.m_leg + p.l^2*p.m_leg + p.l^2*p.mh - p.l^2*p.m_leg*...
         cos(x(1) - x(2))^2));
  dfdu_full_abs(4,1)=-(p.b_len - p.l*cos(x(1) - x(2)))/(p.b_len*(p.a_len^2*p.m_leg + p.l^2*p.m_leg +...
          p.l^2*p.mh - p.l^2*p.m_leg*cos(x(1) - x(2))^2));

 