import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.persistence_sediment_irfmn.ismPersistenceSedimentIrfmn;
import model.ModelExecutionTest;

public class PersistenceSedimentIRFMNTest extends ModelExecutionTest {
    @Override
    protected InsilicoModel getModel() throws InitFailureException, GenericFailureException {
        return new ismPersistenceSedimentIrfmn();
    }
}
