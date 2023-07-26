function PackingFraction(X,ParticleArea)
    return (size(X,2)*ParticleArea) ./ Alphashapes.AlphaShapeVolume(AlphaShape(X,FindAlpha(X)))
end
