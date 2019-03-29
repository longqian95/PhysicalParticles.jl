@inline distance(p1::AbstractPoint, p2::AbstractPoint) = norm(p1-p2)
@inline distance(p1::AbstractParticle, p2::AbstractParticle) = norm(p1.Pos-p2.Pos)
