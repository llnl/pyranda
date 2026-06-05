# Define the equations of motion
euler_1d = """
# Primary Equations of motion here
ddt(:rho:)  =  -ddx(:rho:*:u:)
ddt(:rhou:) =  -ddx(:rhou:*:u: + :p: - :tau:)
ddt(:Et:)   =  -ddx( (:Et: + :p: - :tau:)*:u: )
# Conservative filter of the EoM
:rho:       =  fbar( :rho:  )
:rhou:      =  fbar( :rhou: )
:Et:        =  fbar( :Et:   )
# Update the primatives and enforce the EOS
:u:         =  :rhou: / :rho:
:p:         =  ( :Et: - .5*:rho:*(:u:*:u:) ) * ( :gamma: - 1.0 )
# Artificial bulk viscosity (old school way)
:div:       =  ddx(:u:) 
:beta:      =  gbar(abs(ring(:div:))) * :rho: * 7.0e-2
:tau:       =  :beta:*:div:
"""

euler_2d = """
# Primary Equations of motion here
ddt(:rho:)  =  -ddx(:rho:*:u:)                  - ddy(:rho:*:v:)
ddt(:rhou:) =  -ddx(:rhou:*:u: + :p: - :tau:)   - ddy(:rhou:*:v:)
ddt(:rhov:) =  -ddx(:rhov:*:u:)                 - ddy(:rhov:*:v: + :p: - :tau:)
ddt(:Et:)   =  -ddx( (:Et: + :p: - :tau:)*:u: ) - ddy( (:Et: + :p: - :tau:)*:v: )
# Conservative filter of the EoM
:rho:       =  fbar( :rho:  )
:rhou:      =  fbar( :rhou: )
:rhov:      =  fbar( :rhov: )
:Et:        =  fbar( :Et:   )
# Update the primatives and enforce the EOS
:u:         =  :rhou: / :rho:
:v:         =  :rhov: / :rho:
:p:         =  ( :Et: - .5*:rho:*(:u:*:u: + :v:*:v:) ) * ( :gamma: - 1.0 )
# Artificial bulk viscosity (old school way)
:div:       =  ddx(:u:) + ddy(:v:)
:beta:      =  gbar(abs(ring(:div:))) * :rho: * 7.0e-2
:tau:       =  :beta:*:div:
"""


euler_3d = """
# Primary Equations of motion here
ddt(:rho:)  =  -div(:rho:*:u:  ,  :rho:*:v: , :rho:*:w:)
ddt(:rhou:) =  -div(:rhou:*:u: - :tauxx: , :rhou:*:v: - :tauxy:, :rhou:*:w: - :tauxz:)
ddt(:rhov:) =  -div(:rhov:*:u: - :tauxy: , :rhov:*:v: - :tauyy:, :rhov:*:w: - :tauyz:)
ddt(:rhow:) =  -div(:rhow:*:u: - :tauxz: , :rhow:*:v: - :tauyz:, :rhow:*:w: - :tauzz:)
ddt(:Et:)   =  -div( (:Et: - :tauxx:)*:u: - :tauxy:*:v: - :tauxz:*:w: , (:Et: - :tauyy:)*:v: -:tauxy:*:u: - :tauyz:*:w: ,  (:Et: - :tauzz:)*:w: - :tauxz:*:u: - :tauyz:*:v: )
# Conservative filter of the EoM
:rho:       =  fbar( :rho:  )
:rhou:      =  fbar( :rhou: )
:rhov:      =  fbar( :rhov: )
:rhow:      =  fbar( :rhow: )
:Et:        =  fbar( :Et:   )
# Update the primatives and enforce the EOS
:u:         =  :rhou: / :rho:
:v:         =  :rhov: / :rho:
:w:         =  :rhow: / :rho:
:p:         =  ( :Et: - .5*:rho:*(:u:*:u: + :v:*:v: + :w:*:w:) ) * ( :gamma: - 1.0 )
# Divergence and cross derivatives
:div:       =  div(:u:,:v:,:w:) 
[:ux:,:uy:,:uz:] = grad(:u:)
[:vx:,:vy:,:vz:] = grad(:v:)
[:wx:,:wy:,:wz:] = grad(:w:)
# Artificial bulk viscosity 
:mu:        =  gbar( ringV(:u:,:v:,:w:) * :rho: ) * 1.0e-1 + mu0
:beta:      =  gbar( abs(ring(:div:)) * :rho: )  * 7.0e-2
:taudia:    =  (:beta:-2./3.*:mu:) *:div: - :p:
:tauxx:     =  2.0*:mu:*:ux:   + :taudia:
:tauyy:     =  2.0*:mu:*:vy:   + :taudia:
:tauzz:     =  2.0*:mu:*:wz:   + :taudia:
:tauxy:     = :mu:*(:uy:+:vx:)
:tauxz:     = :mu:*(:uz:+:wx:)
:tauyz:     = :mu:*(:vz:+:wz:)
:cs:  = sqrt( :p: / :rho: * :gamma: )
:dt: = dt.courant(:u:,:v:,:w:,:cs:)*1.0
:dt: = numpy.minimum(:dt:,0.1 * dt.diff(:beta:,:rho:))
:dt: = numpy.minimum(:dt:,0.1 * dt.diff(:mu:,:rho:))
"""


euler_3d_dir = """
# Primary Equations of motion here
ddt(:rho:)  =  -div(:rho:*:u:  ,  :rho:*:v: , :rho:*:w:)
ddt(:rhou:) =  ( -ddx( :FAu:*:detJ: ) - ddy( :FBu:*:detJ: ) - ddz( :FCu:*:detJ: ) ) / :detJ:
ddt(:rhov:) =  ( -ddx( :FAv:*:detJ: ) - ddy( :FBv:*:detJ: ) - ddz( :FCv:*:detJ: ) ) / :detJ:
ddt(:rhow:) =  ( -ddx( :FAw:*:detJ: ) - ddy( :FBw:*:detJ: ) - ddz( :FCw:*:detJ: ) ) / :detJ:
ddt(:Et:)   =  ( -ddx( :FAe:*:detJ: ) - ddy( :FBe:*:detJ: ) - ddz( :FCe:*:detJ: ) ) / :detJ:
# Conservative filter of the EoM
:rho:       =  fbar( :rho:  )
:rhou:      =  fbar( :rhou: )
:rhov:      =  fbar( :rhov: )
:rhow:      =  fbar( :rhow: )
:Et:        =  fbar( :Et:   )
# Update the primatives and enforce the EOS
:u:         =  :rhou: / :rho:
:v:         =  :rhov: / :rho:
:w:         =  :rhow: / :rho:
:p:         =  ( :Et: - .5*:rho:*(:u:*:u: + :v:*:v: + :w:*:w:) ) * ( :gamma: - 1.0 )
# Divergence and cross derivatives
:div:       =  div(:u:,:v:,:w:) 
[:ux:,:uy:,:uz:] = grad(:u:)
[:vx:,:vy:,:vz:] = grad(:v:)
[:wx:,:wy:,:wz:] = grad(:w:)
# Artificial bulk viscosity 
:mu:        =  gbar( ringV(:u:,:v:,:w:) * :rho: ) * 1.0e-4 + mu0
#:beta:      =  gbar( abs(ring(:div:)) * :rho: )  * 7.0e-2
#:mu:        = mu0
:beta:      =  0.0*:beta: #gbar( abs(ring(:div:)) * :rho: )  * 7.0e-2
:taudia:    =  (:beta:-2./3.*:mu:) *:div: - :p:
:tauxx:     =  2.0*:mu:*:ux:   + :taudia:
:tauyy:     =  2.0*:mu:*:vy:   + :taudia:
:tauzz:     =  2.0*:mu:*:wz:   + :taudia:
:tauxy:     = :mu:*(:uy:+:vx:)
:tauxz:     = :mu:*(:uz:+:wx:)
:tauyz:     = :mu:*(:vz:+:wz:)
:rhoA: = ddx(:rho:)
:rhoB: = ddy(:rho:)
:rhoC: = ddz(:rho:)
:rhoM: = sqrt( :rhoA:**2 + :rhoB:**2 + :rhoC:**2 )
:LsA: = abs(:rhoA:)*:dA:/:rhoM:
:LsB: = abs(:rhoB:)*:dB:/:rhoM:
:LsC: = abs(:rhoC:)*:dC:/:rhoM:
:betaT: = gbar(abs(dd4x(:div:)))*:LsA: + gbar(abs(dd4y(:div:)))*:LsB: + gbar(abs(dd4z(:div:)))*:LsC:
:betaA: = :rho: * :betaT: * :dA:
:betaB: = :rho: * :betaT: * :dB:
:betaC: = :rho: * :betaT: * :dC:
:cs:  = sqrt( :p: / :rho: * :gamma: )
:dtC: = dt.courant(:u:,:v:,:w:,:cs:)*1.0
:dt: = numpy.minimum(:dtC:,0.2 * dt.diffDir(:betaA:,:rho:,:dA:))
:dt: = numpy.minimum(:dt: ,0.2 * dt.diffDir(:betaB:,:rho:,:dB:))
:dt: = numpy.minimum(:dt: ,0.2 * dt.diffDir(:betaC:,:rho:,:dC:))
:dt: = numpy.minimum(:dt:,0.1 * dt.diff(:mu:,:rho:))
# Physical fluxes
:Fxu: = :rhou:*:u: - :tauxx:
:Fyu: = :rhou:*:v: - :tauxy:
:Fzu: = :rhou:*:w: - :tauxz:
:Fxv: = :rhov:*:u: - :tauxy:
:Fyv: = :rhov:*:v: - :tauyy:
:Fzv: = :rhov:*:w: - :tauyz:
:Fxw: = :rhow:*:u: - :tauxz:
:Fyw: = :rhow:*:v: - :tauyz:
:Fzw: = :rhow:*:w: - :tauzz:
:Fxe: = (:Et: - :tauxx:)*:u: - :tauxy:*:v: - :tauxz:*:w:
:Fye: = (:Et: - :tauyy:)*:v: - :tauxy:*:u: - :tauyz:*:w:
:Fze: = (:Et: - :tauzz:)*:w: - :tauxz:*:u: - :tauyz:*:v:
# Transformed fluxes + direectional art. bulk terms
:FAu: = :Fxu:*:dAdx: + :Fyu:*:dAdy: + :Fzu:*:dAdz: - :betaA:*:div:*:dAdx:*:detJ:
:FBu: = :Fxu:*:dBdx: + :Fyu:*:dBdy: + :Fzu:*:dBdz: - :betaB:*:div:*:dBdx:*:detJ:
:FCu: = :Fxu:*:dCdx: + :Fyu:*:dCdy: + :Fzu:*:dCdz: - :betaC:*:div:*:dCdx:*:detJ:
:FAv: = :Fxv:*:dAdx: + :Fyv:*:dAdy: + :Fzv:*:dAdz: - :betaA:*:div:*:dAdy:*:detJ:
:FBv: = :Fxv:*:dBdx: + :Fyv:*:dBdy: + :Fzv:*:dBdz: - :betaB:*:div:*:dBdy:*:detJ:
:FCv: = :Fxv:*:dCdx: + :Fyv:*:dCdy: + :Fzv:*:dCdz: - :betaC:*:div:*:dCdy:*:detJ:
:FAw: = :Fxw:*:dAdx: + :Fyw:*:dAdy: + :Fzw:*:dAdz: - :betaA:*:div:*:dAdz:*:detJ:
:FBw: = :Fxw:*:dBdx: + :Fyw:*:dBdy: + :Fzw:*:dBdz: - :betaB:*:div:*:dBdz:*:detJ:
:FCw: = :Fxw:*:dCdx: + :Fyw:*:dCdy: + :Fzw:*:dCdz: - :betaC:*:div:*:dCdz:*:detJ:
:FAe: = :Fxe:*:dAdx: + :Fye:*:dAdy: + :Fze:*:dAdz: - :betaA:*:div:*:detJ:*(:u:*:dAdx:+:v:*:dAdy:+:w:*:dAdz:)
:FBe: = :Fxe:*:dBdx: + :Fye:*:dBdy: + :Fze:*:dBdz: - :betaB:*:div:*:detJ:*(:u:*:dBdx:+:v:*:dBdy:+:w:*:dBdz:)
:FCe: = :Fxe:*:dCdx: + :Fye:*:dCdy: + :Fze:*:dCdz: - :betaC:*:div:*:detJ:*(:u:*:dCdx:+:v:*:dCdy:+:w:*:dCdz:)
"""
