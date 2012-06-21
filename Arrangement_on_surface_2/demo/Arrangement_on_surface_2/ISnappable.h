#ifndef SNAPPABLE_H
#define SNAPPABLE_H

class ISnappable
{
public:
    virtual ~ISnappable( ) { }
    virtual void setSnappingEnabled( bool b ) = 0;
    virtual void setSnapToGridEnabled( bool b ) = 0;
}; // class ISnappable

#endif // SNAPPABLE_H
